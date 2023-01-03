/*
 * Copyright (C) 2010 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

/* This is a JNI example where we use native methods to play sounds
 * using OpenSL ES. See the corresponding Java source file located at:
 *
 *   src/com/example/nativeaudio/NativeAudio/NativeAudio.java
 */
#include <time.h>
#include <android/log.h>
#include <stdlib.h>
#include <assert.h>
#include <jni.h>
#include <string.h>
#include <pthread.h>
#include <fftw3.h>


// for __android_log_print(ANDROID_LOG_INFO, "YourApp", "formatted message");
// #include <android/log.h>

// for native audio
#include <SLES/OpenSLES.h>
#include <SLES/OpenSLES_Android.h>

// for native asset manager
#include <sys/types.h>
#include <android/asset_manager.h>
#include <android/asset_manager_jni.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

// pre-recorded sound clips, both are 8 kHz mono 16-bit signed little endian
static const char hello[] =
#include "hello_clip.h"
;

static const char android[] =
#include "android_clip.h"
;
static long now_ms(void) {
    struct timespec res;
    clock_gettime(CLOCK_REALTIME, &res);
    return (long)(1000.0 * res.tv_sec + (double) res.tv_nsec / 1e6);
}

// engine interfaces
static SLObjectItf engineObject = NULL;
static SLEngineItf engineEngine;

// output mix interfaces
static SLObjectItf outputMixObject = NULL;
static SLEnvironmentalReverbItf outputMixEnvironmentalReverb = NULL;

//int bias=200;
//int win_size=2400;
int filt_offset=127;
int xcorr_counter=0;
int naiser_index_to_process=-1;
int speed=-1;
int sendDelay=-1;
jboolean responder=JNI_FALSE;
int* xcorr_helper(void* context, short* data, int globalOffset, short* globalData);
int* xcorr_helper2(void* context, short* data, int globalOffset, int N);
jboolean freed = JNI_FALSE;

int chirpsPlayed = 0;
double dtx = .8/100;
double drx = 3.3/100;
typedef struct mycontext{
    int totalSegments;
    int queuedSegments;
    int processedSegments;
    int playOffset;
    int bufferSize;
    int timingOffset;
    int initialOffset;
    short* data;
    short* bigdata;
    short* refData;
    char* topfilename;
    char* bottomfilename;
    char* meta_filename;
    double* naiserTx1;
    double* naiserTx2;
    int naiserTx1Count;
    int naiserTx2Count;
    int speed;
    float minPeakDistance;
    jboolean getOneMoreFlag;
    int calibSigLen;
    int sendDelay;
    jboolean runXcorr;
    jboolean sendReply;
    jboolean responder;
    jboolean naiser;
    jint warmdownTime;
    jint preamble_len;
    jfloat xcorrthresh;
    jint N0;
    jboolean CP;
    jfloat naiserThresh;
    jfloat naiserShoulder;
    jint win_size;
    jint bias;
    int seekback;
    double pthresh;
    int recorder_offset;
    int filenum;
    int dataSize;
    float initialDelay;
    long* mic_ts;
    long* speaker_ts;
    int ts_len;
    char* speaker_ts_fname;
    char* mic_ts_fname;
    int bigBufferSize;
    int bigBufferTimes;
}mycontext;

jboolean wroteToDisk=JNI_FALSE;
mycontext* cxt=NULL;
mycontext* cxt2=NULL;
mycontext* cxt3=NULL;
// buffer queue player interfaces
static SLObjectItf bqPlayerObject = NULL;
static SLPlayItf bqPlayerPlay;
static SLAndroidSimpleBufferQueueItf bqPlayerBufferQueue;
static SLEffectSendItf bqPlayerEffectSend;
static SLMuteSoloItf bqPlayerMuteSolo;
static SLVolumeItf bqPlayerVolume;
static SLmilliHertz bqPlayerSampleRate = 0;
static short *resampleBuf = NULL;
jboolean forcewrite=JNI_FALSE;
static int lock=0;
// a mutext to guard against re-entrance to record & playback
// as well as make recording and playing back to be mutually exclusive
// this is to avoid crash at situations like:
//    recording is in session [not finished]
//    user presses record button and another recording coming in
// The action: when recording/playing back is not finished, ignore the new request
static pthread_mutex_t  audioEngineLock = PTHREAD_MUTEX_INITIALIZER;

// aux effect on the output mix, used by the buffer queue player
static const SLEnvironmentalReverbSettings reverbSettings =
        SL_I3DL2_ENVIRONMENT_PRESET_STONECORRIDOR;

// URI player interfaces
static SLObjectItf uriPlayerObject = NULL;

// file descriptor player interfaces
static SLObjectItf fdPlayerObject = NULL;

// recorder interfaces
static SLObjectItf recorderObject = NULL;
static SLRecordItf recorderRecord;
static SLAndroidSimpleBufferQueueItf recorderBufferQueue;

// synthesized sawtooth clip
#define SAWTOOTH_FRAMES 8000
static short sawtoothBuffer[SAWTOOTH_FRAMES];

// 5 seconds of recorded audio at 16 kHz mono, 16-bit signed little endian
jint FS=0;
static short* recorderBuffer;
static unsigned recorderSize = 0;

double distance1=0;
double distance2=0;
double distance3=0;

int self_chirp_idx=-1;
int last_chirp_idx=-1;
double chirp_indexes[5] = {-1,-1,-1,-1,-1};
int lastidx=0;
double last_xcorr_val=0;
double last_naiser_val=0;

jboolean reply_ready=JNI_FALSE;
// pointer and size of the next player buffer to enqueue, and number of remaining buffers
static short *nextBuffer;
static unsigned nextSize;

int receivedIdx=-1;
int replyIdx1 = -1;
int replyIdx2 = -1;
int replyIdx3 = -1;

char* getString_d(double* data, int N) {
    int maxCharLength = 16;
    int maxLineLength = 1023;
    char* str=calloc(maxCharLength*N,sizeof(char));
    int index = 0;
    for (int i = 0; i < N; i++) {
        index += sprintf(&str[index],"%.5f,",(data[i]));
    }
    return str;
}

char* getString_s(short* data, int N) {
    int maxCharLength = 11;
    int maxLineLength = 1023;
    char* str=calloc(maxCharLength*N,sizeof(char));
    int index = 0;
    for (int i = 0; i < N; i++) {
        index += sprintf(&str[index],"%d,",(int)(data[i]));
    }
    return str;
}

// synthesize a mono sawtooth wave and place it into a buffer (called automatically on load)
__attribute__((constructor)) static void onDlOpen(void)
{
    unsigned i;
    for (i = 0; i < SAWTOOTH_FRAMES; ++i) {
        sawtoothBuffer[i] = 32768 - ((i % 100) * 660);
    }
}

void releaseResampleBuf(void) {
    if( 0 == bqPlayerSampleRate) {
        /*
         * we are not using fast path, so we were not creating buffers, nothing to do
         */
        return;
    }

    free(resampleBuf);
    resampleBuf = NULL;
}

/*
 * Only support up-sampling
 */
short* createResampledBuf(uint32_t idx, uint32_t srcRate, unsigned *size) {
    short  *src = NULL;
    short  *workBuf;
    int    upSampleRate;
    int32_t srcSampleCount = 0;

    if(0 == bqPlayerSampleRate) {
        return NULL;
    }
    if(bqPlayerSampleRate % srcRate) {
        /*
         * simple up-sampling, must be divisible
         */
        return NULL;
    }
    upSampleRate = bqPlayerSampleRate / srcRate;

    switch (idx) {
        case 0:
            return NULL;
        case 1: // HELLO_CLIP
            srcSampleCount = sizeof(hello) >> 1;
            src = (short*)hello;
            break;
        case 2: // ANDROID_CLIP
            srcSampleCount = sizeof(android) >> 1;
            src = (short*) android;
            break;
        case 3: // SAWTOOTH_CLIP
            srcSampleCount = SAWTOOTH_FRAMES;
            src = sawtoothBuffer;
            break;
        case 4: // captured frames
            srcSampleCount = recorderSize / sizeof(short);
            src =  recorderBuffer;
            break;
        default:
            assert(0);
            return NULL;
    }

    resampleBuf = (short*) malloc((srcSampleCount * upSampleRate) << 1);
    if(resampleBuf == NULL) {
        return resampleBuf;
    }
    workBuf = resampleBuf;
    for(int sample=0; sample < srcSampleCount; sample++) {
        for(int dup = 0; dup  < upSampleRate; dup++) {
            *workBuf++ = src[sample];
        }
    }

    *size = (srcSampleCount * upSampleRate) << 1;     // sample format is 16 bit
    return resampleBuf;
}

speakerCounter=0;
void static bqPlayerCallback(SLAndroidSimpleBufferQueueItf bq, void *context)
{
    mycontext* cxt = (mycontext*)context;

    cxt->speaker_ts[speakerCounter++] = now_ms();

    assert(bq == bqPlayerBufferQueue);

    if (reply_ready && replyIdx1 >= 0) {
        __android_log_print(ANDROID_LOG_VERBOSE, "debug",
                            "REPLY READY %d %d",receivedIdx,replyIdx1);
        int counter2=0;

        int start_idx = replyIdx1;
        int counter=0;
        for (int i = start_idx; i < start_idx + cxt->preamble_len; i++) {
            cxt->data[i] = cxt->refData[counter++];
        }
        counter2++;
        reply_ready=JNI_FALSE;
        cxt->sendReply=JNI_TRUE;
    }

    SLresult result;
    if (cxt->queuedSegments<cxt->totalSegments) {
        result = (*bqPlayerBufferQueue)->Enqueue(bqPlayerBufferQueue, cxt->data+cxt->playOffset, cxt->bufferSize*sizeof(short));
        assert(SL_RESULT_SUCCESS == result);

//        if (!cxt->responder) {
//            chirpsPlayed=cxt->playOffset/cxt->sendDelay;
//            __android_log_print(ANDROID_LOG_VERBOSE, "debug","CHIRPS PLAYED %d",chirpsPlayed);
//        }
//        else {
//            if (receivedIdx>0) {
//                chirpsPlayed =
//                        (cxt->playOffset - receivedIdx - (cxt->initialDelay * FS)) / cxt->sendDelay;
//                __android_log_print(ANDROID_LOG_VERBOSE, "debug", "CHIRPS PLAYED %d", chirpsPlayed);
//            }
//        }

        cxt->playOffset += cxt->bufferSize;
        cxt->queuedSegments+=1;
    }
    else {
        __android_log_print(ANDROID_LOG_VERBOSE, "debug2", "write speaker %s %d",cxt->speaker_ts_fname,cxt->ts_len);
        FILE* fp = fopen(cxt->speaker_ts_fname,"w+");
        fp = fopen(cxt->speaker_ts_fname,"w+");
        for (int i = 0; i < cxt->ts_len-1; i++) {
            fprintf(fp,"%d\n",cxt->speaker_ts[i]);
        }
        fclose(fp);
        free(cxt->speaker_ts);
        __android_log_print(ANDROID_LOG_VERBOSE, "debug2", "write speaker done");
    }
}

double** fftcomplexoutnative_double(double* data, jint dataN, jint bigN) {
    fftw_complex *in , *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * bigN);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * bigN);

    for (int i = 0; i < bigN; i++) {
        in[i][0] = 0;
        in[i][1] = 0;
        out[i][0] = 0;
        out[i][1] = 0;
    }

    for (int i = 0; i < dataN; i++) {
        in[i][0] = data[i];
    }

    p = fftw_plan_dft_1d(bigN, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

//    __android_log_print(ANDROID_LOG_VERBOSE, "hello", "%f %f %f",data[0],out[0][0],out[1][0]);

    double** outarray = calloc(2,sizeof(double*));
    outarray[0] = calloc(bigN,sizeof(double));
    outarray[1] = calloc(bigN,sizeof(double));

    for (int i = 0; i < bigN; i++) {
        outarray[0][i] = out[i][0];
        outarray[1][i] = out[i][1];
    }

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    return outarray;
}

void conjnative(double** data, int N) {
    for (int i = 0; i < N; i++) {
        data[1][i] = -data[1][i];
    }
}

double** timesnative(double** c1, double** c2, int N) {
    double* real1 = c1[0];
    double* imag1 = c1[1];
    double* real2 = c2[0];
    double* imag2 = c2[1];

    double** outarray = calloc(2,sizeof(double*));
    outarray[0] = calloc(N,sizeof(double));
    outarray[1] = calloc(N,sizeof(double));

    for (int i = 0; i < N; i++) {
        outarray[0][i] = real1[i]*real2[i]-imag1[i]*imag2[i];
        outarray[1][i] = imag1[i]*real2[i]+real1[i]*imag2[i];
    }
//    char* str1=getString_d(outarray[0],N);
//    char* str2=getString_d(outarray[1],N);
    return outarray;
}

double* ifftnative(double** data, int N) {
    fftw_complex *in , *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    for (int i = 0; i < N; i++) {
        in[i][0] = 0;
        in[i][1] = 0;
        out[i][0] = 0;
        out[i][1] = 0;
    }

    double *realArray = data[0];
    double *imagArray = data[1];
    for (int i = 0; i < N; i++) {
        in[i][0] = realArray[i];
    }
    for (int i = 0; i < N; i++) {
        in[i][1] = imagArray[i];
    }

    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    double* realout = calloc(N,sizeof(double));
    int counter = 0;
    for (int i = 0; i < N; i++) {
//    for (int i = N-1; i >= 0; i--) {
        realout[counter++] = out[i][0]/N;
    }

//    char* str5=getString_d(realout,N/2);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    return realout;
}

int* xcorr(double* filteredData, double* refData, int filtLength, int refLength, int i, int globalOffset, int xcorrthresh, double minPeakDistance, int seekback, double pthresh, jboolean getOneMoreFlag) {

//    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***xcorr samples %f %f %f %f %f",
//                        filteredData[0],filteredData[1],filteredData[2],filteredData[3],filteredData[4]);

//    char* str0=getString_d(filteredData,filtLength);
//    char* str1=getString_d(refData,refLength);

//    clock_t begin = clock();
    double** a=fftcomplexoutnative_double(filteredData, filtLength, filtLength);
//    clock_t end = clock();
//    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
//    __android_log_print(ANDROID_LOG_VERBOSE, "time4", "%.4f", time_spent);

//    begin = clock();
    double** b=fftcomplexoutnative_double(refData, refLength, filtLength);
//    end = clock();
//    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
//    __android_log_print(ANDROID_LOG_VERBOSE, "time5", "%.4f", time_spent);

//    char* str1=getString_d(a[0],filtLength);
//    char* str2=getString_d(a[1],filtLength);

    conjnative(b,filtLength);
    double** multout = timesnative(a, b,filtLength);

//    begin = clock();
    double* corr = ifftnative(multout,filtLength);
//    end = clock();
//    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
//    __android_log_print(ANDROID_LOG_VERBOSE, "time6", "%.4f", time_spent);

    free(a[0]);
    free(a[1]);
    free(b[0]);
    free(b[1]);
    free(multout[0]);
    free(multout[1]);

//    char* str2=getString_d(a[0],filtLength);
//    char* str3=getString_d(b[0],filtLength);
//    char* str4=getString_d(multout[0],filtLength);
//    char* str5=getString_d(corr,filtLength);

    double maxval=-1;
    int maxidx=-1;
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug4","exceed %d",filtLength);
    for (int i = 0; i < filtLength; i++) {
        if (corr[i] > maxval) {
//            __android_log_print(ANDROID_LOG_VERBOSE, "debug4","exceed %f %f",corr[i],maxval);
            maxval=corr[i];
            maxidx=i;
        }
    }
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug4","xcorr3 %f %d %d",
//                        maxval,xcorrthresh,maxidx);
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug","naive xcorr %d %f",maxidx,maxval);

//    if (self_chirp_idx==-1) {
        // seekback, we advance forward in android as the xcorr result is flipped
        for (int i = seekback; i >= 0; i--) {
            int l_new=maxidx-i;
            if (l_new >= 0 && l_new < filtLength) {
                if (corr[l_new] > maxval * pthresh) {
                    maxidx = l_new;
                    break;
                }
            }
        }
//    }

    free(corr);

//    __android_log_print(ANDROID_LOG_VERBOSE, "debug","seekback xcorr %d %f",maxidx,maxval);

//    maxidx = filtLength-(maxidx*2)-1;

//    __android_log_print(ANDROID_LOG_VERBOSE, "hello", "peak %.0f %.0f",maxval,maxval/1000000);
//    maxval /= 1000000;
    last_xcorr_val=maxval;

    // return local index, max val
    int* out = calloc(2,sizeof(int));

    int globalidx=maxidx+globalOffset;
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug","xcorr mindist %d %d %d %d %d %f",maxidx,globalOffset,globalidx,lastidx,
//                        maxidx+globalOffset - lastidx,minPeakDistance);
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug","xcorr mindist %d %d %d %d %d",maxval > xcorrthresh,globalidx - lastidx > minPeakDistance*FS,lastidx==0,globalidx==lastidx,getOneMoreFlag);

    if (maxval > xcorrthresh && (globalidx - lastidx > minPeakDistance*FS||lastidx==0||globalidx==lastidx||getOneMoreFlag)) {
        out[0] = maxidx;
        out[1] = maxval;
//        lastidx=maxidx+globalOffset;
        return out;
    }
    else {
        out[0] = -1;
        out[1] = maxval;
        return out;
    }
}

double* myfir_s(short* data,double* h, int lenData, int lenH) {
    int nconv = lenH+lenData-1;

    double* temp = calloc(nconv,sizeof(double));

    for (int i=0; i<nconv; i++) {
        temp[i]=0;
    }
    for (int n = 0; n < nconv; n++){
        int kmin, kmax;

        temp[n] = 0;

        kmin = (n >= lenH - 1) ? n - (lenH - 1) : 0;
        kmax = (n < lenData - 1) ? n : lenData - 1;

        for (int k = kmin; k <= kmax; k++) {
            temp[n] += data[k] * h[n - k];
        }
    }

    return temp;
}

double* myfir_d(double* data,double* h, int lenData, int lenH) {
    int nconv = lenH+lenData-1;

    double* temp = calloc(nconv,sizeof(double));

    for (int i=0; i<nconv; i++) {
        temp[i]=0;
    }
    for (int n = 0; n < nconv; n++){
        int kmin, kmax;

        temp[n] = 0;

        kmin = (n >= lenH - 1) ? n - (lenH - 1) : 0;
        kmax = (n < lenData - 1) ? n : lenData - 1;

        for (int k = kmin; k <= kmax; k++) {
            temp[n] += data[k] * h[n - k];
        }
    }

    return temp;
}

double sum_multiple(double* seg1, double* seg2, double seg1len, double seg2len) {
    double out = 0;
    for (int i = 0; i < seg1len; i++) {
        out += seg1[i] * seg2[i];
    }
    return out;
}

double* filter_s(short* data, int N) {
//    double h[129]={0.000182981959336989,0.000281596242979397,0.000278146045925432,0.000175593510303356,-4.62343304809110e-19,-0.000200360419689133,-0.000360951066953434,-0.000412991328953827,-0.000300653514269367,1.50969613210241e-18,0.000465978441137468,0.00102258147190892,0.00155366635073282,0.00192784355834049,0.00203610307379658,0.00183119111593167,0.00135603855040258,0.000749280342023432,0.000220958680427305,-2.30998316092680e-19,0.000264751350126403,0.00107567949517942,0.00233229627684471,0.00377262258405226,0.00502311871694995,0.00569227410155340,0.00548603174290836,0.00431270688583350,0.00234309825736876,0,-0.00213082326080853,-0.00345315860003227,-0.00353962427138790,-0.00228740687362881,3.04497082814396e-18,0.00264042059820580,0.00471756289619728,0.00531631859929438,0.00379212348265709,-1.99178361782373e-17,-0.00558925290045386,-0.0119350370250123,-0.0176440086545525,-0.0213191620435004,-0.0219592383672450,-0.0193021014035209,-0.0140079988312124,-0.00761010753938763,-0.00221480214028708,-3.05997847316821e-18,-0.00261988600696407,-0.0106615919949044,-0.0233019683686087,-0.0382766989959210,-0.0522058307502582,-0.0612344107673335,-0.0618649720104803,-0.0518011604210193,-0.0306049219949348,2.05563094310381e-17,0.0362729247397242,0.0730450660051671,0.104645449260725,0.125975754818137,0.133506082359307,0.125975754818137,0.104645449260725,0.0730450660051671,0.0362729247397242,2.05563094310381e-17,-0.0306049219949348,-0.0518011604210193,-0.0618649720104803,-0.0612344107673335,-0.0522058307502582,-0.0382766989959210,-0.0233019683686087,-0.0106615919949044,-0.00261988600696407,-3.05997847316821e-18,-0.00221480214028708,-0.00761010753938763,-0.0140079988312124,-0.0193021014035209,-0.0219592383672450,-0.0213191620435004,-0.0176440086545525,-0.0119350370250123,-0.00558925290045386,-1.99178361782373e-17,0.00379212348265709,0.00531631859929438,0.00471756289619728,0.00264042059820580,3.04497082814396e-18,-0.00228740687362881,-0.00353962427138790,-0.00345315860003227,-0.00213082326080853,0,0.00234309825736876,0.00431270688583350,0.00548603174290836,0.00569227410155340,0.00502311871694995,0.00377262258405226,0.00233229627684471,0.00107567949517942,0.000264751350126403,-2.30998316092680e-19,0.000220958680427305,0.000749280342023432,0.00135603855040258,0.00183119111593167,0.00203610307379658,0.00192784355834049,0.00155366635073282,0.00102258147190892,0.000465978441137468,1.50969613210241e-18,-0.000300653514269367,-0.000412991328953827,-0.000360951066953434,-0.000200360419689133,-4.62343304809110e-19,0.000175593510303356,0.000278146045925432,0.000281596242979397,0.000182981959336989};
    double h[257]={-7.75418371456968e-05,8.07613875877651e-06,-3.92526167043826e-05,-0.000196577687451030,-0.000363815429427690,-0.000429011103355452,-0.000344211593290773,-0.000161433405990234,-2.27906993160878e-06,2.21324905722762e-05,-0.000110576061155958,-0.000311661348347424,-0.000434696883769424,-0.000373945469366597,-0.000142686343390440,0.000125526059860333,0.000260911644546323,0.000175459780613039,-6.58953350088675e-05,-0.000276653523425456,-0.000270803670261422,2.76171636530950e-06,0.000409094778853374,0.000705093329509247,0.000703814715600538,0.000415348528110552,5.95769502689691e-05,-7.20325819790529e-05,0.000178836072582780,0.000707880805548229,0.00119598585719552,0.00131796961455718,0.000977094933481562,0.000397469014108926,-7.25350872688717e-06,8.40942390412198e-05,0.000657376794616464,0.00133985891012699,0.00163482563945913,0.00126878481436442,0.000410120695812505,-0.000411562869857602,-0.000655022323947711,-0.000157555976487184,0.000713059390223053,0.00127654946849575,0.00100368894695367,-8.78349944396206e-05,-0.00139727891686581,-0.00212196774639508,-0.00181669359290785,-0.000743094372055135,0.000260732796658281,0.000336747387397122,-0.000775766687272849,-0.00250158378691291,-0.00377229131416849,-0.00375110821922860,-0.00244950041720240,-0.000783946821672791,1.22028724702117e-05,-0.000749781225902183,-0.00267602888721496,-0.00448654067178731,-0.00486180048480386,-0.00338257879880225,-0.000896858215572822,0.000992067265165680,0.00101589067172485,-0.000826917138183567,-0.00317607429010071,-0.00418112513485569,-0.00278805783020865,0.000447875236091432,0.00363839200129231,0.00480975817205362,0.00330192421198659,0.000331793636070950,-0.00174191471604462,-0.00102206334558371,0.00253562740495695,0.00695068161204548,0.00947644401597403,0.00849636178646048,0.00475354493893717,0.000989204169992018,0.000130947535939305,0.00321428726881982,0.00849734359572208,0.0124875619151081,0.0122927532000494,0.00772334699889114,0.00163146724296435,-0.00188914557070325,-0.000327202667360082,0.00532653020483520,0.0110004458963365,0.0121917731549891,0.00714175294801697,-0.00167720293411387,-0.00896015759588886,-0.0101054903654515,-0.00456104382870571,0.00341275770579504,0.00733760762480312,0.00293769160556652,-0.00857852574844589,-0.0208135125902682,-0.0262005495076531,-0.0212830286030305,-0.00964933176381504,-0.000190413898523792,-0.00137902871462172,-0.0151351588915629,-0.0345440205887962,-0.0475222390535327,-0.0446652617466912,-0.0264003884645968,-0.00439713781596734,0.00461270470496120,-0.0103659562820506,-0.0455949628374347,-0.0815323553384488,-0.0915801894325124,-0.0566750715438293,0.0222527355408236,0.121454331780568,0.203704303656161,0.235571062399159,0.203704303656161,0.121454331780568,0.0222527355408236,-0.0566750715438293,-0.0915801894325124,-0.0815323553384488,-0.0455949628374347,-0.0103659562820506,0.00461270470496120,-0.00439713781596734,-0.0264003884645968,-0.0446652617466912,-0.0475222390535327,-0.0345440205887962,-0.0151351588915629,-0.00137902871462172,-0.000190413898523792,-0.00964933176381504,-0.0212830286030305,-0.0262005495076531,-0.0208135125902682,-0.00857852574844589,0.00293769160556652,0.00733760762480312,0.00341275770579504,-0.00456104382870571,-0.0101054903654515,-0.00896015759588886,-0.00167720293411387,0.00714175294801697,0.0121917731549891,0.0110004458963365,0.00532653020483520,-0.000327202667360082,-0.00188914557070325,0.00163146724296435,0.00772334699889114,0.0122927532000494,0.0124875619151081,0.00849734359572208,0.00321428726881982,0.000130947535939305,0.000989204169992018,0.00475354493893717,0.00849636178646048,0.00947644401597403,0.00695068161204548,0.00253562740495695,-0.00102206334558371,-0.00174191471604462,0.000331793636070950,0.00330192421198659,0.00480975817205362,0.00363839200129231,0.000447875236091432,-0.00278805783020865,-0.00418112513485569,-0.00317607429010071,-0.000826917138183567,0.00101589067172485,0.000992067265165680,-0.000896858215572822,-0.00338257879880225,-0.00486180048480386,-0.00448654067178731,-0.00267602888721496,-0.000749781225902183,1.22028724702117e-05,-0.000783946821672791,-0.00244950041720240,-0.00375110821922860,-0.00377229131416849,-0.00250158378691291,-0.000775766687272849,0.000336747387397122,0.000260732796658281,-0.000743094372055135,-0.00181669359290785,-0.00212196774639508,-0.00139727891686581,-8.78349944396206e-05,0.00100368894695367,0.00127654946849575,0.000713059390223053,-0.000157555976487184,-0.000655022323947711,-0.000411562869857602,0.000410120695812505,0.00126878481436442,0.00163482563945913,0.00133985891012699,0.000657376794616464,8.40942390412198e-05,-7.25350872688717e-06,0.000397469014108926,0.000977094933481562,0.00131796961455718,0.00119598585719552,0.000707880805548229,0.000178836072582780,-7.20325819790529e-05,5.95769502689691e-05,0.000415348528110552,0.000703814715600538,0.000705093329509247,0.000409094778853374,2.76171636530950e-06,-0.000270803670261422,-0.000276653523425456,-6.58953350088675e-05,0.000175459780613039,0.000260911644546323,0.000125526059860333,-0.000142686343390440,-0.000373945469366597,-0.000434696883769424,-0.000311661348347424,-0.000110576061155958,2.21324905722762e-05,-2.27906993160878e-06,-0.000161433405990234,-0.000344211593290773,-0.000429011103355452,-0.000363815429427690,-0.000196577687451030,-3.92526167043826e-05,8.07613875877651e-06,-7.75418371456968e-05};
    double* filtered = myfir_s(data,h,N,257);
    return filtered;
}

double* filter_d(double* data, int N) {
//    double h[129]={0.000182981959336989,0.000281596242979397,0.000278146045925432,0.000175593510303356,-4.62343304809110e-19,-0.000200360419689133,-0.000360951066953434,-0.000412991328953827,-0.000300653514269367,1.50969613210241e-18,0.000465978441137468,0.00102258147190892,0.00155366635073282,0.00192784355834049,0.00203610307379658,0.00183119111593167,0.00135603855040258,0.000749280342023432,0.000220958680427305,-2.30998316092680e-19,0.000264751350126403,0.00107567949517942,0.00233229627684471,0.00377262258405226,0.00502311871694995,0.00569227410155340,0.00548603174290836,0.00431270688583350,0.00234309825736876,0,-0.00213082326080853,-0.00345315860003227,-0.00353962427138790,-0.00228740687362881,3.04497082814396e-18,0.00264042059820580,0.00471756289619728,0.00531631859929438,0.00379212348265709,-1.99178361782373e-17,-0.00558925290045386,-0.0119350370250123,-0.0176440086545525,-0.0213191620435004,-0.0219592383672450,-0.0193021014035209,-0.0140079988312124,-0.00761010753938763,-0.00221480214028708,-3.05997847316821e-18,-0.00261988600696407,-0.0106615919949044,-0.0233019683686087,-0.0382766989959210,-0.0522058307502582,-0.0612344107673335,-0.0618649720104803,-0.0518011604210193,-0.0306049219949348,2.05563094310381e-17,0.0362729247397242,0.0730450660051671,0.104645449260725,0.125975754818137,0.133506082359307,0.125975754818137,0.104645449260725,0.0730450660051671,0.0362729247397242,2.05563094310381e-17,-0.0306049219949348,-0.0518011604210193,-0.0618649720104803,-0.0612344107673335,-0.0522058307502582,-0.0382766989959210,-0.0233019683686087,-0.0106615919949044,-0.00261988600696407,-3.05997847316821e-18,-0.00221480214028708,-0.00761010753938763,-0.0140079988312124,-0.0193021014035209,-0.0219592383672450,-0.0213191620435004,-0.0176440086545525,-0.0119350370250123,-0.00558925290045386,-1.99178361782373e-17,0.00379212348265709,0.00531631859929438,0.00471756289619728,0.00264042059820580,3.04497082814396e-18,-0.00228740687362881,-0.00353962427138790,-0.00345315860003227,-0.00213082326080853,0,0.00234309825736876,0.00431270688583350,0.00548603174290836,0.00569227410155340,0.00502311871694995,0.00377262258405226,0.00233229627684471,0.00107567949517942,0.000264751350126403,-2.30998316092680e-19,0.000220958680427305,0.000749280342023432,0.00135603855040258,0.00183119111593167,0.00203610307379658,0.00192784355834049,0.00155366635073282,0.00102258147190892,0.000465978441137468,1.50969613210241e-18,-0.000300653514269367,-0.000412991328953827,-0.000360951066953434,-0.000200360419689133,-4.62343304809110e-19,0.000175593510303356,0.000278146045925432,0.000281596242979397,0.000182981959336989};
    double h[257]={-7.75418371456968e-05,8.07613875877651e-06,-3.92526167043826e-05,-0.000196577687451030,-0.000363815429427690,-0.000429011103355452,-0.000344211593290773,-0.000161433405990234,-2.27906993160878e-06,2.21324905722762e-05,-0.000110576061155958,-0.000311661348347424,-0.000434696883769424,-0.000373945469366597,-0.000142686343390440,0.000125526059860333,0.000260911644546323,0.000175459780613039,-6.58953350088675e-05,-0.000276653523425456,-0.000270803670261422,2.76171636530950e-06,0.000409094778853374,0.000705093329509247,0.000703814715600538,0.000415348528110552,5.95769502689691e-05,-7.20325819790529e-05,0.000178836072582780,0.000707880805548229,0.00119598585719552,0.00131796961455718,0.000977094933481562,0.000397469014108926,-7.25350872688717e-06,8.40942390412198e-05,0.000657376794616464,0.00133985891012699,0.00163482563945913,0.00126878481436442,0.000410120695812505,-0.000411562869857602,-0.000655022323947711,-0.000157555976487184,0.000713059390223053,0.00127654946849575,0.00100368894695367,-8.78349944396206e-05,-0.00139727891686581,-0.00212196774639508,-0.00181669359290785,-0.000743094372055135,0.000260732796658281,0.000336747387397122,-0.000775766687272849,-0.00250158378691291,-0.00377229131416849,-0.00375110821922860,-0.00244950041720240,-0.000783946821672791,1.22028724702117e-05,-0.000749781225902183,-0.00267602888721496,-0.00448654067178731,-0.00486180048480386,-0.00338257879880225,-0.000896858215572822,0.000992067265165680,0.00101589067172485,-0.000826917138183567,-0.00317607429010071,-0.00418112513485569,-0.00278805783020865,0.000447875236091432,0.00363839200129231,0.00480975817205362,0.00330192421198659,0.000331793636070950,-0.00174191471604462,-0.00102206334558371,0.00253562740495695,0.00695068161204548,0.00947644401597403,0.00849636178646048,0.00475354493893717,0.000989204169992018,0.000130947535939305,0.00321428726881982,0.00849734359572208,0.0124875619151081,0.0122927532000494,0.00772334699889114,0.00163146724296435,-0.00188914557070325,-0.000327202667360082,0.00532653020483520,0.0110004458963365,0.0121917731549891,0.00714175294801697,-0.00167720293411387,-0.00896015759588886,-0.0101054903654515,-0.00456104382870571,0.00341275770579504,0.00733760762480312,0.00293769160556652,-0.00857852574844589,-0.0208135125902682,-0.0262005495076531,-0.0212830286030305,-0.00964933176381504,-0.000190413898523792,-0.00137902871462172,-0.0151351588915629,-0.0345440205887962,-0.0475222390535327,-0.0446652617466912,-0.0264003884645968,-0.00439713781596734,0.00461270470496120,-0.0103659562820506,-0.0455949628374347,-0.0815323553384488,-0.0915801894325124,-0.0566750715438293,0.0222527355408236,0.121454331780568,0.203704303656161,0.235571062399159,0.203704303656161,0.121454331780568,0.0222527355408236,-0.0566750715438293,-0.0915801894325124,-0.0815323553384488,-0.0455949628374347,-0.0103659562820506,0.00461270470496120,-0.00439713781596734,-0.0264003884645968,-0.0446652617466912,-0.0475222390535327,-0.0345440205887962,-0.0151351588915629,-0.00137902871462172,-0.000190413898523792,-0.00964933176381504,-0.0212830286030305,-0.0262005495076531,-0.0208135125902682,-0.00857852574844589,0.00293769160556652,0.00733760762480312,0.00341275770579504,-0.00456104382870571,-0.0101054903654515,-0.00896015759588886,-0.00167720293411387,0.00714175294801697,0.0121917731549891,0.0110004458963365,0.00532653020483520,-0.000327202667360082,-0.00188914557070325,0.00163146724296435,0.00772334699889114,0.0122927532000494,0.0124875619151081,0.00849734359572208,0.00321428726881982,0.000130947535939305,0.000989204169992018,0.00475354493893717,0.00849636178646048,0.00947644401597403,0.00695068161204548,0.00253562740495695,-0.00102206334558371,-0.00174191471604462,0.000331793636070950,0.00330192421198659,0.00480975817205362,0.00363839200129231,0.000447875236091432,-0.00278805783020865,-0.00418112513485569,-0.00317607429010071,-0.000826917138183567,0.00101589067172485,0.000992067265165680,-0.000896858215572822,-0.00338257879880225,-0.00486180048480386,-0.00448654067178731,-0.00267602888721496,-0.000749781225902183,1.22028724702117e-05,-0.000783946821672791,-0.00244950041720240,-0.00375110821922860,-0.00377229131416849,-0.00250158378691291,-0.000775766687272849,0.000336747387397122,0.000260732796658281,-0.000743094372055135,-0.00181669359290785,-0.00212196774639508,-0.00139727891686581,-8.78349944396206e-05,0.00100368894695367,0.00127654946849575,0.000713059390223053,-0.000157555976487184,-0.000655022323947711,-0.000411562869857602,0.000410120695812505,0.00126878481436442,0.00163482563945913,0.00133985891012699,0.000657376794616464,8.40942390412198e-05,-7.25350872688717e-06,0.000397469014108926,0.000977094933481562,0.00131796961455718,0.00119598585719552,0.000707880805548229,0.000178836072582780,-7.20325819790529e-05,5.95769502689691e-05,0.000415348528110552,0.000703814715600538,0.000705093329509247,0.000409094778853374,2.76171636530950e-06,-0.000270803670261422,-0.000276653523425456,-6.58953350088675e-05,0.000175459780613039,0.000260911644546323,0.000125526059860333,-0.000142686343390440,-0.000373945469366597,-0.000434696883769424,-0.000311661348347424,-0.000110576061155958,2.21324905722762e-05,-2.27906993160878e-06,-0.000161433405990234,-0.000344211593290773,-0.000429011103355452,-0.000363815429427690,-0.000196577687451030,-3.92526167043826e-05,8.07613875877651e-06,-7.75418371456968e-05};
    double* filtered = myfir_d(data,h,N,257);
    return filtered;
}

void setReply(int idx, mycontext* cxt) {
    // don't trigger on the calibration chirp
//    __android_log_print(ANDROID_LOG_VERBOSE,"debug2","set reply start %d %d",idx,replyIdx1);
//    if (cxt->sendReply&&cxt->responder && idx > FS && idx-replyIdx1 > FS) {
    if (cxt->responder && idx > FS && idx-replyIdx1 > FS) {
//    if (cxt->responder) {
        receivedIdx = idx;
//        replyIdx1 = idx + cxt->initialDelay*FS;
        replyIdx1 = idx - cxt->timingOffset + cxt->sendDelay;
//        replyIdx2 = replyIdx1 + cxt->sendDelay;
//        replyIdx3 = replyIdx2 + cxt->sendDelay;
//        if (replyIdx3 < cxt->dataSize-cxt->preamble_len-cxt->bufferSize &&
//            self_chirp_idx != -1 && idx != self_chirp_idx) {
        if (replyIdx1 < cxt->dataSize-cxt->preamble_len-cxt->bufferSize) {
//            __android_log_print(ANDROID_LOG_VERBOSE, "debug", "send index %d %d %d %d %f",
//                                idx, cxt->timingOffset, cxt->sendDelay, replyIdx1,
//                                (replyIdx1 - idx) / (double) FS);
//            __android_log_print(ANDROID_LOG_VERBOSE,"debug2","set reply %d %d %d",
//                                cxt->initialDelay, idx, replyIdx1);
            reply_ready = JNI_TRUE;
            cxt->sendReply = JNI_FALSE;
            chirpsPlayed+=1;
        }
    }
}

jboolean timeOffsetUpdated=JNI_FALSE;
void updateTimingOffset(int global_xcorr_idx, int local_chirp_idx, mycontext* cxt) {
//    __android_log_print(ANDROID_LOG_VERBOSE,"debug2","update timing offset");
    double oneSampleDelay = 1.0/FS;
    int transmitDelay = (int) (((33.0 / 1000) / cxt->speed)/oneSampleDelay);
    cxt->timingOffset = global_xcorr_idx - (cxt->initialOffset) - transmitDelay;
    if (cxt->timingOffset < 0) {
//        __android_log_print(ANDROID_LOG_VERBOSE, "debug","break");
    }
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug",
//                        "UPDATE TIMING OFFSET %d %d %d %d %d %d %d",global_xcorr_idx, local_chirp_idx, cxt->initialOffset,
//                        transmitDelay,cxt->timingOffset,cxt->speed,cxt->processedSegments);
    timeOffsetUpdated=JNI_TRUE;
}

void* xcorr_thread(void* context) {
    mycontext* cxt = (mycontext*)context;
//    __android_log_print(ANDROID_LOG_VERBOSE, "hello", "thread created");
    if (cxt->runXcorr) {
        int global_xcorr_idx=-1;
        int local_xcorr_idx=-1;

        if (!cxt->getOneMoreFlag) {
//            __android_log_print(ANDROID_LOG_VERBOSE, "debug",
//                                "xcorr helper1 %d %d %d",cxt->processedSegments,(cxt->processedSegments) * cxt->bigBufferSize, cxt->bigBufferSize);
            int *result = NULL;
            int globalOffset=0;
            if (cxt->processedSegments == 0) {
                short *data = cxt->data + (cxt->bigBufferSize * (cxt->processedSegments));
//                char* str1=getString_d(data,cxt->bigBufferSize);
//                __android_log_print(ANDROID_LOG_VERBOSE, "debug5","offset1 %d",(cxt->bigBufferSize * (cxt->processedSegments)));
                globalOffset = (cxt->processedSegments) * cxt->bigBufferSize;
                result = xcorr_helper2(context, data, globalOffset, cxt->bigBufferSize);
            }
            else {
                short *data = cxt->data + (cxt->bigBufferSize * (cxt->processedSegments-1));
//                char* str1=getString_d(data,cxt->bigBufferSize);
//                __android_log_print(ANDROID_LOG_VERBOSE, "debug5","offset2 %d",(cxt->bigBufferSize * (cxt->processedSegments-1)));
                globalOffset = (cxt->processedSegments-1) * cxt->bigBufferSize;
                result = xcorr_helper2(context, data, globalOffset, cxt->bigBufferSize*2);
            }

            local_xcorr_idx = result[0];
            int naiser_idx = result[2];

            if (local_xcorr_idx>=0) {
                int synclag = cxt->seekback;
//                __android_log_print(ANDROID_LOG_VERBOSE, "debug", "get flag? %d %d %d %d",
//                                    local_xcorr_idx, local_xcorr_idx + cxt->calibSigLen + synclag,
//                                    local_xcorr_idx + cxt->naiserTx2Count + cxt->win_size,
//                                    cxt->bigBufferSize);

                jboolean c1 = cxt->processedSegments == 0 && local_xcorr_idx + cxt->naiserTx2Count + synclag >= cxt->bigBufferSize;
                jboolean c2 = cxt->processedSegments == 0 && local_xcorr_idx + cxt->naiserTx2Count + cxt->win_size >= cxt->bigBufferSize;
                jboolean c3 = cxt->processedSegments > 0 && local_xcorr_idx + cxt->naiserTx2Count + synclag >= cxt->bigBufferSize*2;
                jboolean c4 = cxt->processedSegments > 0 && local_xcorr_idx + cxt->naiserTx2Count + cxt->win_size >= cxt->bigBufferSize*2;
                if (c1 || c2 || c3 || c4) {
//                    __android_log_print(ANDROID_LOG_VERBOSE, "debug",
//                                        "get one more flag set to true");
                    cxt->getOneMoreFlag = JNI_TRUE;
                } else if (naiser_idx >= 0){
                    global_xcorr_idx = result[0] + globalOffset;
                    if (xcorr_counter<5) {
//                        __android_log_print(ANDROID_LOG_VERBOSE, "chirp", "got chirp1 %d %d %d",xcorr_counter,global_xcorr_idx,naiser_idx);
                        chirp_indexes[xcorr_counter++] = global_xcorr_idx;
                        lastidx=global_xcorr_idx;
                    }
                    last_chirp_idx = global_xcorr_idx;
//                    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "xcorr1 %d %d %d %d %d %d",
//                                        local_xcorr_idx,
//                                        global_xcorr_idx, result[1], result[2], globalOffset,
//                                        xcorr_counter);

//                    if (cxt->sendReply&&cxt->responder&&xcorr_counter==2) {
                    if (cxt->responder) {
//                        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "set reply2");
                        setReply(global_xcorr_idx, cxt);
                    }
                }
                else {
//                    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "get flag is false");
                }
            }

            free(result);
        }
        else {
            // look back half a second
            int globalOffset = (cxt->processedSegments-1) * (cxt->bigBufferSize);

            short *data = cxt->data +
                          (int) ((cxt->bigBufferSize * (double) (cxt->processedSegments - 1)));

//            __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***samples %d %hi %hi %hi %hi %hi",
//                                cxt->bigBufferSize * (cxt->processedSegments),data[0],data[1],data[2],data[3],data[4]);

//            __android_log_print(ANDROID_LOG_VERBOSE, "debug",
//                                "xcorr helper2 %d %d %d",cxt->processedSegments,(cxt->processedSegments-1) * cxt->bigBufferSize,cxt->bigBufferSize);

            int *result = xcorr_helper2(context, data, globalOffset, cxt->bigBufferSize*2);

            local_xcorr_idx = result[0];
            int naiser_out = result[2];

            if (local_xcorr_idx >= 0 && naiser_out >= 0) {
                global_xcorr_idx = result[0]+globalOffset;
                if (xcorr_counter<5) {
//                    __android_log_print(ANDROID_LOG_VERBOSE, "chirp", "got chirp2 %d %d %d",xcorr_counter,global_xcorr_idx,naiser_out);
                    chirp_indexes[xcorr_counter++] = global_xcorr_idx;
                    lastidx=global_xcorr_idx;
                }
                last_chirp_idx = global_xcorr_idx;
//                __android_log_print(ANDROID_LOG_VERBOSE, "debug", "xcorr2 %d %d %d %d", local_xcorr_idx,
//                                    global_xcorr_idx, result[1], xcorr_counter);

//                __android_log_print(ANDROID_LOG_VERBOSE, "debug", "get one more flag set to false 1");
                cxt->getOneMoreFlag = JNI_FALSE;

//                if (cxt->sendReply&&cxt->responder&&xcorr_counter==2) {
                if (cxt->responder) {
//                    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "set reply");
                    setReply(global_xcorr_idx, cxt);
                }
            } else if (local_xcorr_idx < 0 || naiser_out < 0) {
                // occurs in the case of noise
//                __android_log_print(ANDROID_LOG_VERBOSE, "debug", "get one more flag set to false 2");
                cxt->getOneMoreFlag = JNI_FALSE;
            }
            free(result);
        }

        // assume that self chirp is always within first second
        // force an update if it's xcorr or naiser refinement
        if (!cxt->getOneMoreFlag &&
            (cxt->processedSegments == 0 ||
            (global_xcorr_idx<FS && local_xcorr_idx >= 0 && global_xcorr_idx >= 0))) {
            self_chirp_idx = global_xcorr_idx;
            updateTimingOffset(global_xcorr_idx,local_xcorr_idx,cxt);
        }
    }
    cxt->processedSegments+=1;
}

double getdist(int earier_chirp_idx, int later_chirp_index,int delay,int dtx,int drx, double speed) {
    double diff = (later_chirp_index - earier_chirp_idx) - delay;
    double delta = diff / (double)FS * speed;
    double distance = (delta / 2) + dtx + drx;
    return distance;
}

void getdistall() {
    if (!responder) {
        distance1 = -1;
        distance2 = -1;
        distance3 = -1;

        if (chirp_indexes[1]!=-1) {
            distance1 = getdist(self_chirp_idx, chirp_indexes[1], sendDelay, dtx, drx,
                                speed);
        }
        if (chirp_indexes[2]!=-1) {
            distance2 = getdist(self_chirp_idx, chirp_indexes[2], sendDelay * 2, dtx, drx,
                                speed);
        }
        if (chirp_indexes[3]!=-1) {
            distance3 = getdist(self_chirp_idx, chirp_indexes[3], sendDelay * 3, dtx,
                                drx, speed);
        }
    }
}

JNIEXPORT void JNICALL
Java_com_example_nativeaudio_NativeAudio_forcewrite(JNIEnv *env, jclass clazz) {
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "forcewrite");
    forcewrite=JNI_TRUE;
    getdistall();
}

void stopithelper() {
    if (!freed) {
//        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "stopithelper");
        SLresult result;
        result = (*recorderRecord)->SetRecordState(recorderRecord, SL_RECORDSTATE_STOPPED);
        assert(SL_RESULT_SUCCESS == result);

        result = (*bqPlayerPlay)->SetPlayState(bqPlayerPlay, SL_PLAYSTATE_STOPPED);
        assert(SL_RESULT_SUCCESS == result);

//        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "freeing1");
//    free(chirp_indexes);
        if (cxt != NULL) {
            free(cxt);
            free(cxt->data);
            free(cxt->refData);
        }
//        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "freeing2");
        if (cxt2 != NULL) {
            free(cxt2);
            free(cxt2->data);
        }
//        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "freeing3");
        if (cxt3 != NULL) {
            free(cxt3);
            free(cxt3->data);
        }
//        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "done freeing");
        freed = JNI_TRUE;
    }
}

int smallBufferIdx=0;
int micCounter=0;
// this callback handler is called every time a buffer finishes recording
void bqRecorderCallback(SLAndroidSimpleBufferQueueItf bq, void *context)
{
    assert(bq == recorderBufferQueue);

    mycontext* cxt=(mycontext*)context;
    cxt->mic_ts[micCounter++] = now_ms();

    int offset=0;
    // determine whether top/bottom microphone is the odd/even sample in the buffer
    if (smallBufferIdx==0) {
        int maxval1=0;
        int maxval2=0;

        int counter=0;
        for (int i = 0; i < cxt->bufferSize; i++) {
            int val1 = cxt->bigdata[(cxt->processedSegments*cxt->bufferSize*2)+counter];
            int val2 = cxt->bigdata[(cxt->processedSegments*cxt->bufferSize*2)+counter+1];
            if (val1 > maxval1) {
                maxval1=val1;
            }
            if (val2 > maxval2) {
                maxval2=val2;
            }
            counter+=2;
        }
        if (maxval1 > maxval2) {
            cxt->recorder_offset = 0;
        }
        else {
            cxt->recorder_offset = 1;
        }
    }

    // we only process data from one microphone, we use the one with higher amplitude (probably bottom one)
    int counter=cxt->recorder_offset;
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug5", "copy %d %d %d %d",smallBufferIdx,cxt->bufferSize,smallBufferIdx*cxt->bufferSize,smallBufferIdx*cxt->bufferSize*2);
    for (int i = 0; i < cxt->bufferSize; i++) {
        cxt->data[smallBufferIdx*cxt->bufferSize+i] = cxt->bigdata[(smallBufferIdx*cxt->bufferSize*2)+counter];
        counter+=2;
    }
    smallBufferIdx+=1;

    if (cxt->queuedSegments<cxt->totalSegments) {
        SLresult result=(*recorderBufferQueue)->Enqueue(recorderBufferQueue,
                                                        (cxt->bigdata)+(cxt->bufferSize*2*cxt->queuedSegments),
                                                        cxt->bufferSize*2*sizeof(short));
        assert(SL_RESULT_SUCCESS == result);
        cxt->queuedSegments+=1;
    }

    // bug: if the time is too short, then this will use 'old' values of an index (i.e. not naiser-refined)
    jboolean c1 = forcewrite&&!wroteToDisk;
    jboolean c2 = !cxt->responder && xcorr_counter==4;
    jboolean c3 = !cxt->responder&&cxt->queuedSegments==cxt->totalSegments;
    jboolean c5 = (cxt->responder&&cxt->queuedSegments==cxt->totalSegments);
    if ((c1 || c2 || c3 || c5) && !wroteToDisk) {
//        __android_log_print(ANDROID_LOG_VERBOSE, "debug2", "write to mic disk %s %d",cxt->mic_ts_fname,cxt->ts_len);

        // write to disk
        wroteToDisk=JNI_TRUE;
        FILE* fp = fopen(cxt->mic_ts_fname,"w+");
        fp = fopen(cxt->mic_ts_fname,"w+");
        for (int i = 0; i < cxt->ts_len-1; i++) {
            fprintf(fp,"%d\n",cxt->mic_ts[i]);
        }
        fclose(fp);
        free(cxt->mic_ts);

        fp = fopen(cxt->bottomfilename,"w+");
        for (int i = cxt->recorder_offset; i < cxt->totalSegments*cxt->bufferSize*2-2; i+=2) {
            fprintf(fp,"%d ",cxt->bigdata[i]);
        }
        fclose(fp);

        fp = fopen(cxt->topfilename,"w+");
        for (int i = -(cxt->recorder_offset-1); i < cxt->totalSegments*cxt->bufferSize*2-2; i+=2) {
            fprintf(fp,"%d ",cxt->bigdata[i]);
        }
        fclose(fp);

        ///////////////////////////////////////////////////////////////////////

        wroteToDisk=JNI_TRUE;

//        __android_log_print(ANDROID_LOG_VERBOSE, "debug2", "finish writing");
        stopithelper();
    }

    if (!wroteToDisk && cxt->runXcorr) {
        if (micCounter%cxt->bigBufferTimes==0) {
//            __android_log_print(ANDROID_LOG_VERBOSE, "debug2", "process %d %d %d",micCounter,cxt->bigBufferTimes,
//                                cxt->mic_ts[micCounter-1]);
//            char* str1=getString_s(cxt->bigdata,cxt->bigBufferSize);
//            char* str2=getString_s(cxt->data,cxt->bigBufferSize);
            pthread_t t;
            pthread_create(&t, NULL, xcorr_thread, context);
        }
    }
}

// create the engine and output mix objects
JNIEXPORT jdoubleArray JNICALL
Java_com_example_nativeaudio_NativeAudio_getDistance(JNIEnv* env, jclass clazz, jboolean reply)
{
    if (reply) {
        jdouble out[5];
        out[0] = self_chirp_idx;
        out[1] = receivedIdx;
        out[2] = replyIdx1;
        out[3] = replyIdx2;
        out[4] = replyIdx3;

        jdoubleArray result;
        result = (*env)->NewDoubleArray(env, 5);
        (*env)->SetDoubleArrayRegion(env, result, 0, 5, out);

        return result;
    }
    else {
        jdouble out[7];
        out[0] = distance1;
        out[1] = distance2;
        out[2] = distance3;
        out[3] = chirp_indexes[0];
        out[4] = chirp_indexes[1];
        out[5] = chirp_indexes[2];
        out[6] = chirp_indexes[3];

        jdoubleArray result;
        result = (*env)->NewDoubleArray(env, 7);
        (*env)->SetDoubleArrayRegion(env, result, 0, 7, out);

        return result;
    }
}

// create the engine and output mix objects
JNIEXPORT jdoubleArray JNICALL
Java_com_example_nativeaudio_NativeAudio_getVal(JNIEnv* env, jclass clazz)
{
    jdouble out[2];
    out[0] = last_xcorr_val;
    out[1] = last_naiser_val;

    jdoubleArray result;
    result = (*env)->NewDoubleArray(env,2);
    (*env)->SetDoubleArrayRegion(env, result, 0, 2, out);

    return result;
}
// create the engine and output mix objects
JNIEXPORT jboolean JNICALL
Java_com_example_nativeaudio_NativeAudio_responderDone(JNIEnv* env, jclass clazz)
{
    return cxt2->responder&&!cxt2->sendReply&&wroteToDisk;
}

JNIEXPORT jboolean JNICALL
Java_com_example_nativeaudio_NativeAudio_replySet(JNIEnv* env, jclass clazz)
{
    return !cxt2->sendReply;
}

JNIEXPORT jint JNICALL
Java_com_example_nativeaudio_NativeAudio_getXcorrCount(JNIEnv* env, jclass clazz)
{
    return xcorr_counter;
}

JNIEXPORT jintArray JNICALL
Java_com_example_nativeaudio_NativeAudio_getReplyIndexes(JNIEnv* env, jclass clazz)
{
    jint out[3];
    out[0] = replyIdx1;
    out[1] = replyIdx2;
    out[2] = replyIdx3;

    jintArray result;
    result = (*env)->NewIntArray(env,3);
    (*env)->SetIntArrayRegion(env, result, 0, 3, out);

    return result;
}

JNIEXPORT jint JNICALL
Java_com_example_nativeaudio_NativeAudio_getQueuedSpeakerSegments(JNIEnv* env, jclass clazz)
{
    return cxt->queuedSegments;
}

// create the engine and output mix objects
JNIEXPORT void JNICALL
Java_com_example_nativeaudio_NativeAudio_createEngine(JNIEnv* env, jclass clazz)
{
    SLresult result;

    // create engine
    result = slCreateEngine(&engineObject, 0, NULL, 0, NULL, NULL);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;

    // realize the engine
    result = (*engineObject)->Realize(engineObject, SL_BOOLEAN_FALSE);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;

    // get the engine interface, which is needed in order to create other objects
    result = (*engineObject)->GetInterface(engineObject, SL_IID_ENGINE, &engineEngine);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;

    // create output mix, with environmental reverb specified as a non-required interface
    const SLInterfaceID ids[1] = {SL_IID_ENVIRONMENTALREVERB};
    const SLboolean req[1] = {SL_BOOLEAN_FALSE};
    result = (*engineEngine)->CreateOutputMix(engineEngine, &outputMixObject, 1, ids, req);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;

    // realize the output mix
    result = (*outputMixObject)->Realize(outputMixObject, SL_BOOLEAN_FALSE);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;

    // get the environmental reverb interface
    // this could fail if the environmental reverb effect is not available,
    // either because the feature is not present, excessive CPU load, or
    // the required MODIFY_AUDIO_SETTINGS permission was not requested and granted
    result = (*outputMixObject)->GetInterface(outputMixObject, SL_IID_ENVIRONMENTALREVERB,
                                              &outputMixEnvironmentalReverb);
    if (SL_RESULT_SUCCESS == result) {
        result = (*outputMixEnvironmentalReverb)->SetEnvironmentalReverbProperties(
                outputMixEnvironmentalReverb, &reverbSettings);
        (void)result;
    }
    // ignore unsuccessful result codes for environmental reverb, as it is optional for this example
}

// create buffer queue audio player
JNIEXPORT void JNICALL
Java_com_example_nativeaudio_NativeAudio_createBufferQueueAudioPlayer(JNIEnv* env, jclass clazz, jint sampleRate, jint bufSize)
{
    SLresult result;
    if (sampleRate >= 0 && bufSize >= 0 ) {
        bqPlayerSampleRate = sampleRate * 1000;
        /*
         * device native buffer size is another factor to minimize audio latency, not used in this
         * sample: we only play one giant buffer here
         */
    }

    // configure audio source
    int numBuffers=5;
    SLDataLocator_AndroidSimpleBufferQueue loc_bufq = {SL_DATALOCATOR_ANDROIDSIMPLEBUFFERQUEUE, numBuffers};
    SLDataFormat_PCM format_pcm = {SL_DATAFORMAT_PCM, 1, SL_SAMPLINGRATE_8,
                                   SL_PCMSAMPLEFORMAT_FIXED_16, SL_PCMSAMPLEFORMAT_FIXED_16,
                                   SL_SPEAKER_FRONT_CENTER, SL_BYTEORDER_LITTLEENDIAN};
    /*
     * Enable Fast Audio when possible:  once we set the same rate to be the native, fast audio path
     * will be triggered
     */
    if(bqPlayerSampleRate) {
        format_pcm.samplesPerSec = bqPlayerSampleRate;       //sample rate in mili second
    }
    SLDataSource audioSrc = {&loc_bufq, &format_pcm};

    // configure audio sink
    SLDataLocator_OutputMix loc_outmix = {SL_DATALOCATOR_OUTPUTMIX, outputMixObject};
    SLDataSink audioSnk = {&loc_outmix, NULL};

    /*
     * create audio player:
     *     fast audio does not support when SL_IID_EFFECTSEND is required, skip it
     *     for fast audio case
     */
    const SLInterfaceID ids[3] = {SL_IID_BUFFERQUEUE, SL_IID_VOLUME,
//                                  SL_IID_EFFECTSEND,
            /*SL_IID_MUTESOLO,*/};
    const SLboolean req[3] = {SL_BOOLEAN_TRUE, SL_BOOLEAN_TRUE, SL_BOOLEAN_TRUE,
            /*SL_BOOLEAN_TRUE,*/ };

    result = (*engineEngine)->CreateAudioPlayer(engineEngine, &bqPlayerObject, &audioSrc, &audioSnk,
                                                bqPlayerSampleRate? 2 : 3, ids, req);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;

    // realize the player
    result = (*bqPlayerObject)->Realize(bqPlayerObject, SL_BOOLEAN_FALSE);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;

    // get the play interface
    result = (*bqPlayerObject)->GetInterface(bqPlayerObject, SL_IID_PLAY, &bqPlayerPlay);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;

    // get the buffer queue interface
    result = (*bqPlayerObject)->GetInterface(bqPlayerObject, SL_IID_BUFFERQUEUE,
                                             &bqPlayerBufferQueue);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;

    // get the effect send interface
    bqPlayerEffectSend = NULL;
    if( 0 == bqPlayerSampleRate) {
        result = (*bqPlayerObject)->GetInterface(bqPlayerObject, SL_IID_EFFECTSEND,
                                                 &bqPlayerEffectSend);
        assert(SL_RESULT_SUCCESS == result);
        (void)result;
    }

#if 0   // mute/solo is not supported for sources that are known to be mono, as this is
    // get the mute/solo interface
    result = (*bqPlayerObject)->GetInterface(bqPlayerObject, SL_IID_MUTESOLO, &bqPlayerMuteSolo);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;
#endif

    // get the volume interface
    result = (*bqPlayerObject)->GetInterface(bqPlayerObject, SL_IID_VOLUME, &bqPlayerVolume);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;
}

// create audio recorder: recorder is not in fast path
//    like to avoid excessive re-sampling while playing back from Hello & Android clip
JNIEXPORT jboolean JNICALL
Java_com_example_nativeaudio_NativeAudio_createAudioRecorder(JNIEnv* env, jclass clazz)
{
    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "start recorder setup");

    SLresult result;

    int numchannels=2;
    int mics=0;
    if (numchannels==1) {
        mics=SL_SPEAKER_FRONT_CENTER;
    }
    else {
        mics=SL_SPEAKER_FRONT_LEFT | SL_SPEAKER_FRONT_RIGHT;
    }

    // configure audio source
    SLDataLocator_IODevice loc_dev = {SL_DATALOCATOR_IODEVICE, SL_IODEVICE_AUDIOINPUT,
                                      SL_DEFAULTDEVICEID_AUDIOINPUT, NULL};
    SLDataSource audioSrc = {&loc_dev, NULL};

    SLDataLocator_AndroidSimpleBufferQueue loc_bq = {SL_DATALOCATOR_ANDROIDSIMPLEBUFFERQUEUE, 1};
    SLDataFormat_PCM format_pcm = {SL_DATAFORMAT_PCM, numchannels, SL_SAMPLINGRATE_44_1,
                                   SL_PCMSAMPLEFORMAT_FIXED_16, SL_PCMSAMPLEFORMAT_FIXED_16,
                                   mics, SL_BYTEORDER_LITTLEENDIAN};
    if (FS==48000) {
        format_pcm.samplesPerSec=SL_SAMPLINGRATE_48;
    }

    SLDataSink audioSnk = {&loc_bq, &format_pcm};

    // create audio recorder
    // (requires the RECORD_AUDIO permission)
    const SLInterfaceID id[2] = {SL_IID_ANDROIDCONFIGURATION, SL_IID_ANDROIDSIMPLEBUFFERQUEUE};
    const SLboolean req[1] = {SL_BOOLEAN_TRUE};
    result = (*engineEngine)->CreateAudioRecorder(engineEngine, &recorderObject, &audioSrc,
                                                  &audioSnk, 2, id, req);
    if (SL_RESULT_SUCCESS != result) {
        return JNI_FALSE;
    }

    // Configure the voice recognition preset which has no
    // signal processing for lower latency.
    SLAndroidConfigurationItf inputConfig;
    result = (*recorderObject)
            ->GetInterface(recorderObject, SL_IID_ANDROIDCONFIGURATION,
                           &inputConfig);

    if (SL_RESULT_SUCCESS == result) {
        SLuint32 presetValue = SL_ANDROID_RECORDING_PRESET_CAMCORDER;
//        SLuint32 presetValue = SL_ANDROID_RECORDING_PRESET_VOICE_RECOGNITION;
        (*inputConfig)
                ->SetConfiguration(inputConfig, SL_ANDROID_KEY_RECORDING_PRESET,
                                   &presetValue, sizeof(SLuint32));
    }

    // realize the audio recorder
    result = (*recorderObject)->Realize(recorderObject, SL_BOOLEAN_FALSE);
    if (SL_RESULT_SUCCESS != result) {
        return JNI_FALSE;
    }

    // get the record interface
    result = (*recorderObject)->GetInterface(recorderObject, SL_IID_RECORD, &recorderRecord);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;

    // get the buffer queue interface
    result = (*recorderObject)->GetInterface(recorderObject, SL_IID_ANDROIDSIMPLEBUFFERQUEUE,
                                             &recorderBufferQueue);
    assert(SL_RESULT_SUCCESS == result);
    (void)result;

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "finish recorder setup");

    return JNI_TRUE;
}

void shutdownhelper() {
    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "shutdown1");
    // destroy buffer queue audio player object, and invalidate all associated interfaces
    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "destroy player1");

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "destroy player2");
    // destroy file descriptor audio player object, and invalidate all associated interfaces
    if (fdPlayerObject != NULL) {
        (*fdPlayerObject)->Destroy(fdPlayerObject);
        fdPlayerObject = NULL;
    }

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "destroy player3");
    // destroy URI audio player object, and invalidate all associated interfaces
    if (uriPlayerObject != NULL) {
        (*uriPlayerObject)->Destroy(uriPlayerObject);
        uriPlayerObject = NULL;
    }

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "destroy recorder");

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "destroy outputmix");
    // destroy output mix object, and invalidate all associated interfaces
    if (outputMixObject != NULL) {
        (*outputMixObject)->Destroy(outputMixObject);
        outputMixObject = NULL;
        outputMixEnvironmentalReverb = NULL;
    }

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "destroy engine");
    // destroy engine object, and invalidate all associated interfaces
    if (lock==1) {
        pthread_mutex_destroy(&audioEngineLock);
    }

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "shutdown2");
}

// shut down the native audio system
JNIEXPORT void JNICALL
Java_com_example_nativeaudio_NativeAudio_shutdown(JNIEnv* env, jclass clazz)
{
    shutdownhelper();
}

JNIEXPORT void JNICALL
Java_com_example_nativeaudio_NativeAudio_stopit(JNIEnv *env, jclass clazz) {
    stopithelper();
}

JNIEXPORT void JNICALL
Java_com_example_nativeaudio_NativeAudio_reset(JNIEnv *env, jclass clazz) {
    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "reset");

    SLresult result;
    result=(*bqPlayerPlay)->SetPlayState(bqPlayerPlay, SL_PLAYSTATE_STOPPED);
    assert(SL_RESULT_SUCCESS == result);

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "reset 1");
    result=(*bqPlayerBufferQueue)->Clear(bqPlayerBufferQueue);
    assert(SL_RESULT_SUCCESS == result);
    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "reset 2");

//    if (recorderRecord!=NULL) {
        result = (*recorderRecord)->SetRecordState(recorderRecord, SL_RECORDSTATE_STOPPED);
//        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "reset rec %d",result);
        assert(SL_RESULT_SUCCESS == result);
//    }

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "reset 3");
    result = (*recorderBufferQueue)->Clear(recorderBufferQueue);
    assert(SL_RESULT_SUCCESS == result);

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "reset done");
}
/////////////////////////////////////////////////////////////////////////////////////////

JNIEXPORT void JNICALL
Java_com_example_nativeaudio_NativeAudio_calibrate(JNIEnv *env, jclass clazz,jshortArray tempData,
                                                   jshortArray tempRefData, jint bSize, jint recordTime,
                                                   jstring ttfilename, jstring tbfilename, jstring tmeta_filename,
                                                   jint initialOffset, jint warmdownTime,
                                                   jint preamble_len, jboolean water,
                                                   jboolean reply, jboolean naiser,
                                                   jint tempSendDelay, jfloat xcorrthresh, jfloat minPeakDistance,
                                                   jint fs, jdoubleArray tnaiserTx1, jdoubleArray tnaiserTx2,
                                                   jint N0, jboolean CP, jfloat naiserThresh, jfloat naiserShoulder,
                                                   jint win_size, jint bias, jint seekback, jdouble pthresh, int round,
                                                   int filenum, jboolean runxcorr, jfloat initialDelay, jstring mic_ts_fname,
                                                   jstring speaker_ts_fname, int bigBufferSize,int bigBufferTimes) {
    freed=JNI_FALSE;
    timeOffsetUpdated=JNI_FALSE;
    int round0 = 0;
    smallBufferIdx=0;
    micCounter=0;
    speakerCounter=0;
    FS=fs;
    chirpsPlayed=0;
    forcewrite=JNI_FALSE;
    responder = reply;
    int bufferSize = bSize;
    reply_ready=JNI_FALSE;
    receivedIdx=-1;
    replyIdx1=-1;
    replyIdx2=-1;
    replyIdx3=-1;
    SLresult result;
    recorderSize = 0;
    self_chirp_idx=-1;
    last_chirp_idx=-1;
    naiser_index_to_process=-1;

    jint N = (*env)->GetArrayLength(env, tempData);
    jint N_ref = (*env)->GetArrayLength(env, tempRefData);

    jint naiserTx1Count = (*env)->GetArrayLength(env, tnaiserTx1);
    jint naiserTx2Count = (*env)->GetArrayLength(env, tnaiserTx2);
    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "naiser %d %d", naiserTx1Count,naiserTx2Count);

    double* naiserTx1 = (double*)((*env)->GetDoubleArrayElements(env, tnaiserTx1, NULL));
    double* naiserTx2 = (double*)((*env)->GetDoubleArrayElements(env, tnaiserTx2, NULL));

    int totalSpeakerLoops = (recordTime*FS)/bufferSize;
    int totalRecorderLoops = (recordTime*FS)/bufferSize;

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "total speaker/recorder loop %d %d", totalSpeakerLoops,totalRecorderLoops);

    short* data = (short*)((*env)->GetShortArrayElements(env, tempData, NULL));
    short* refData = (short*)((*env)->GetShortArrayElements(env, tempRefData, NULL));

    if (round==round0) {
        cxt = calloc(1, sizeof(mycontext));
        cxt->sendDelay=tempSendDelay;
        cxt->preamble_len = preamble_len;
        cxt->data=calloc(bufferSize * totalSpeakerLoops, sizeof(short));
        cxt->initialDelay = initialDelay;
        cxt->refData=calloc(N_ref, sizeof(short));
        memcpy(cxt->refData,refData, N_ref*sizeof(short));

        // copy at beginning
        memcpy(cxt->data,data, N*sizeof(short));

        // copy N seconds later
//        if (!runxcorr && responder) {
//            replyIdx1 = initialDelay*fs;
//            replyIdx2 = replyIdx1 + tempSendDelay;
//            replyIdx3 = replyIdx2 + tempSendDelay;
//
//            int counter = 0;
//            for (int i = replyIdx1; i < replyIdx1+N; i++) {
//                cxt->data[i] = cxt->data[counter++];
//            }
//            counter=0;
//            for (int i = replyIdx2; i < replyIdx2+N; i++) {
//                cxt->data[i] = cxt->data[counter++];
//            }
//            counter=0;
//            for (int i = replyIdx3; i < replyIdx3+N; i++) {
//                cxt->data[i] = cxt->data[counter++];
//            }
//        }
        if (!runxcorr && !responder){
            for (int index = 0; index < bufferSize * totalSpeakerLoops; index += tempSendDelay) {
                __android_log_print(ANDROID_LOG_VERBOSE, "debug", "FILLING %d %d",index, bufferSize * totalSpeakerLoops);
                for (int i = 0; i < N_ref; i++) {
                    cxt->data[index+i] = refData[i];
                }
            }
        }
    }
    else {
        memset(cxt->data,0,bufferSize * totalSpeakerLoops*sizeof(short));
    }

    cxt->warmdownTime=warmdownTime;
    cxt->initialOffset=initialOffset;
    cxt->bufferSize=bufferSize;
    cxt->playOffset=bufferSize;
    cxt->totalSegments=totalSpeakerLoops;
    cxt->queuedSegments=1;
    cxt->responder=reply;
    char* speaker_ts_filename_str = (*env)->GetStringUTFChars(env, speaker_ts_fname, NULL);
    cxt->speaker_ts_fname=speaker_ts_filename_str;
    cxt->speaker_ts = calloc((recordTime*FS)/bufferSize,sizeof(long));
    cxt->ts_len=(recordTime*FS)/bufferSize;

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "populate cxt");

    if (round==round0) {
        result = (*bqPlayerBufferQueue)->RegisterCallback(bqPlayerBufferQueue, bqPlayerCallback, cxt);
        assert(SL_RESULT_SUCCESS == result); //this fails
        (void)result;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    char* topfilename = (*env)->GetStringUTFChars(env, ttfilename, NULL);
    char* bottomfilename = (*env)->GetStringUTFChars(env, tbfilename, NULL);
    char* meta_filename = (*env)->GetStringUTFChars(env, tmeta_filename, NULL);
    char* mic_ts_filename_str = (*env)->GetStringUTFChars(env, mic_ts_fname, NULL);

    wroteToDisk=JNI_FALSE;
    lastidx=0;
    last_xcorr_val=0;
    last_naiser_val=0;
    distance1=-1;
    distance2=-1;
    distance3=-1;
    xcorr_counter=0;
    //cxt is for the speaker
    //cxt2 is for the microphone
    if (round==round0) {
        cxt2 = calloc(1, sizeof(mycontext));
        cxt2->bigdata = calloc(bufferSize * 2 * totalRecorderLoops, sizeof(short));
        cxt2->data = calloc(bufferSize * totalRecorderLoops, sizeof(short));
        cxt2->dataSize = bufferSize*totalRecorderLoops;
    }
    else {
        memset(cxt2->bigdata,0,bufferSize * 2 * totalRecorderLoops * sizeof(short));
        memset(cxt2->data,0,bufferSize * totalRecorderLoops * sizeof(short));
    }
    cxt2->ts_len=(recordTime*FS)/bufferSize;
    cxt2->bigBufferSize=bigBufferSize;
    cxt2->bigBufferTimes=bigBufferTimes;
    cxt2->mic_ts = calloc((recordTime*FS)/bufferSize,sizeof(long));
    cxt2->initialDelay = initialDelay;
    cxt2->xcorrthresh=xcorrthresh;
    cxt2->naiserTx1=naiserTx1;
    cxt2->naiserTx2=naiserTx2;
    cxt2->filenum=filenum;
    cxt2->naiserTx1Count=naiserTx1Count;
    cxt2->naiserTx2Count=naiserTx2Count;
    cxt2->refData=refData;
    cxt2->recorder_offset=0;
    cxt2->totalSegments=totalRecorderLoops;
    cxt2->naiser=naiser;
    cxt2->preamble_len = preamble_len;
    cxt2->topfilename = topfilename;
    cxt2->processedSegments=0;
    cxt2->bottomfilename = bottomfilename;
    cxt2->meta_filename = meta_filename;
    cxt2->mic_ts_fname=mic_ts_filename_str;
    cxt2->bufferSize=bufferSize;
    cxt2->queuedSegments=1;
    cxt2->initialOffset=initialOffset;
    cxt2->warmdownTime = warmdownTime;
    cxt2->naiserThresh=naiserThresh;
    cxt2->naiserShoulder=naiserShoulder;
    cxt2->N0=N0;
    cxt2->CP=CP;
    cxt2->win_size=win_size;
    cxt2->bias=bias;
    cxt2->seekback=seekback;
    cxt2->pthresh=pthresh;

    if (water) {
        cxt2->speed = 1500;
    }
    else {
        cxt2->speed=340;
    }
    cxt2->getOneMoreFlag = JNI_FALSE;
    cxt2->sendDelay=tempSendDelay;
    cxt2->minPeakDistance=minPeakDistance;
    cxt2->runXcorr=runxcorr;
    cxt2->sendReply=JNI_TRUE;
    cxt2->responder=reply;
    speed=cxt2->speed;
    sendDelay=cxt2->sendDelay;
    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "populate cxt2");

    for (int i = 0; i < 5; i++) {
        chirp_indexes[i]=-1;
    }

    if (round==round0) {
        cxt3 = calloc(1, sizeof(mycontext));
        cxt3->data = calloc(bufferSize * totalRecorderLoops, sizeof(short));
    }
    else {
        memset(cxt3->data,0,bufferSize * totalRecorderLoops*sizeof(short));
    }
    cxt3->bufferSize=bufferSize;
    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "populate cxt3");

    if (round==round0) {
        result = (*recorderBufferQueue)->RegisterCallback(recorderBufferQueue, bqRecorderCallback,
                                                          cxt2);
        assert(SL_RESULT_SUCCESS == result);
    }

    result = (*recorderBufferQueue)->Enqueue(recorderBufferQueue, cxt2->bigdata, bufferSize * 2 * sizeof(short));
    assert(SL_RESULT_SUCCESS == result);

    result = (*bqPlayerBufferQueue)->Enqueue(bqPlayerBufferQueue, cxt->data, bufferSize * sizeof(short));
    assert(SL_RESULT_SUCCESS == result);

    result=(*recorderRecord)->SetRecordState(recorderRecord, SL_RECORDSTATE_RECORDING);
    assert(SL_RESULT_SUCCESS == result);

//    usleep(30*1000);

    result = (*bqPlayerPlay)->SetPlayState(bqPlayerPlay, SL_PLAYSTATE_PLAYING);
    assert(SL_RESULT_SUCCESS == result);

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "finish setup");
}

JNIEXPORT void JNICALL
Java_com_example_nativeaudio_NativeAudio_testFileWrite(JNIEnv *env, jclass clazz, jstring tfilename) {
    char* filename = (*env)->GetStringUTFChars(env, tfilename, NULL);
    FILE* fp = fopen(filename,"w+");
    clock_t begin = clock();

    for (int i = 0; i < 48000*30; i++) {
        fprintf(fp,"%d, ",i);
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "time %f", time_spent);

    fclose(fp);
}

// z = abs(c)
double Complex_abs(fftw_complex c){
    return sqrt(pow(c[0], 2 ) + pow(c[1], 2 ) );
}

// z = power(c, 2)
double Complex_power(fftw_complex c){
    return pow(c[0], 2 ) + pow(c[1], 2 );
}

void abs_complex(double * out, fftw_complex* in, unsigned long int Nu){
    unsigned long int  i = 0;
    for( i = 0; i < Nu ; ++i){
        out[i] = Complex_abs(in[i]);
    }
}

void estimate_H(fftw_complex* H, fftw_complex* X, fftw_complex * Y, int fft_begin, int fft_end, unsigned long int Nu){
    unsigned long int  i = 0;
    for( i = 0; i < Nu ; ++i){
        if(i >= fft_begin && i < fft_end){
            double y_pow = Complex_power(X[i]);
            H[i][0] = (Y[i][0]*X[i][0] + X[i][1]*Y[i][1])/y_pow;
            H[i][1] = (Y[i][1]*X[i][0] - Y[i][0]*X[i][1])/y_pow;
        }
        else{
            H[i][0] = 0;
            H[i][1] = 0;
        }
    }
}

void channel_estimation_freq_single(fftw_complex* H_result, double * tx, double * rx, unsigned long int sig_len, int BW1, int BW2, int fs){

    fftw_complex *rx_fft=(fftw_complex*)fftw_malloc(sig_len*sizeof(fftw_complex));
    if(rx_fft==NULL){fprintf(stderr,"malloc failed\n");exit(1);}

    fftw_complex *tx_fft=(fftw_complex*)fftw_malloc(sig_len*sizeof(fftw_complex));
    if(tx_fft==NULL){fprintf(stderr,"malloc failed\n");exit(1);}


    fftw_plan p1=fftw_plan_dft_r2c_1d(sig_len, tx, tx_fft, FFTW_ESTIMATE);
    if (p1==NULL){fprintf(stderr,"plan creation failed\n");exit(1);}

    fftw_plan p2=fftw_plan_dft_r2c_1d(sig_len, rx, rx_fft, FFTW_ESTIMATE);
    if (p2==NULL){fprintf(stderr,"plan creation failed\n");exit(1);}


    fftw_execute(p1);
    fftw_execute(p2);
    //cout <<sig_len <<endl;

    double delta_f = (double)fs/(double)sig_len;
    int fft_begin = (int) round(BW1/delta_f);
    int fft_end = (int) round(BW2/delta_f);

    estimate_H(H_result, tx_fft, rx_fft, fft_begin, fft_end, sig_len);

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);

    fftw_free(tx_fft);
    fftw_free(rx_fft);

    //delete []abs_h;
}


double sum_sqr(double* in, int begin, int end){
    double sum0 = 0;
    for(int i = begin; i < end; ++i){
        sum0 += in[i]*in[i];
    }
    return sum0;
}


double multiptle_sum(double* in1, int begin1, int end1, double * in2 ,int begin2, int end2){
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "msum %d %d %d %d",begin1,end1,begin2,end2);
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "msum %d %d",end1-begin1,end2-begin2);
    assert(end1-begin1 == end2 - begin2);
    double sum0 = 0;
    for(int i = 0; i < end1 - begin1 ; ++i){
        sum0 += in1[begin1 + i]*in2[begin2+i];
    }

    return sum0;
}

int find_max(double* in, unsigned long int len_in){
//    int maxPosition = max_element(in,in+len_in) - in;
    double maxval=0;
    int maxidx=0;
    for (int i = 0; i < len_in; i++) {
        if (in[i] > maxval) {
            maxval=in[i];
            maxidx=i;
        }
    }
//    int maxPosition=maxval-in;
    return maxidx;
}

// received signal, total_length is length of signal
// Nu=length of symbol (720)
// N0 = guard interval (480), we will use CP version
// DIVIDE_FACTOR = 2
// include_zero=false (false for cp version, true for guard interval version)
// threshold=.6
int naiser_corr(double* signal, int total_length , int Nu, int N0, int DIVIDE_FACTOR, jboolean include_zero, double threshold, double nshoulder){
    __android_log_print(ANDROID_LOG_VERBOSE, "debug2", "naiser corr %d",total_length);

    short PN_seq[8] = {1, -1, -1, -1, -1, -1, 1, -1};

    int N_both = Nu + N0;
    int preamble_L = 8*N_both;

    if(total_length-preamble_L < 0){
//        cout << "input signal too short" << endl;
        __android_log_print(ANDROID_LOG_VERBOSE, "debug","input signal too short");
        return -1;
    }

    int len_corr = (total_length - preamble_L)/DIVIDE_FACTOR + 2;
    int num = 0;
//    double* Mn = new double[len_corr];
    double* Mn = calloc(len_corr,sizeof(double));
    memset(Mn, 0, sizeof(double)*len_corr);

    for(int i = 0; i < total_length - preamble_L - 1; i += DIVIDE_FACTOR){
        //cout << i << ' ' << total_length - preamble_L - 1 << ' '  << len_corr  << endl;
//        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "it %d/%d",i,total_length - preamble_L - 1);
        double Pd = 0;
        for(int k = 0; k < 8 - 1; ++k){
//            __android_log_print(ANDROID_LOG_VERBOSE, "debug", "it2 %d",k);
            int bk = PN_seq[k]*PN_seq[k+1];
            double temp_P = multiptle_sum(signal + i, k*N_both + N0, (k + 1)*N_both,
                                          signal + i, (k + 1)*N_both + N0, (k+2)*N_both);
            Pd +=  temp_P*bk;
        }
        double Rd = 0;
//        if(include_zero){
//            Rd = (sum_sqr(signal + i, 0, preamble_L)*Nu)/N_both;
//        }
//        else{
            for(int k = 0; k < 8; ++k){
                Rd += sum_sqr(signal + i, k*N_both + N0, (k + 1)*N_both);
            }
//        }

        Mn[num] = Pd/Rd;

        num ++;
    }

    int max_idx = find_max(Mn, (unsigned long int)len_corr);
    double max_value = Mn[max_idx];

//    char* str=getString_d(Mn,len_corr);
//    char* str2=getString_d(signal,total_length);

    last_naiser_val = max_value;
    __android_log_print(ANDROID_LOG_VERBOSE, "debug6","naiser value %.2f",max_value);
    if(max_value < threshold){
        __android_log_print(ANDROID_LOG_VERBOSE, "debug2","naiser_err1 %.2f %.2f",max_value,threshold);
//        delete [] Mn;
        free(Mn);
        return -1;
    }
    else {
        __android_log_print(ANDROID_LOG_VERBOSE, "debug2","naiser_pass1 %.2f %.2f",max_value,threshold);
    }

    double shoulder = nshoulder*max_value;
    int right = -1;
    int left = -1;

    // find the right shoulder point
    for(int i = max_idx; i < len_corr - 1; ++i){
        if(Mn[i] >= shoulder && Mn[i+1] <= shoulder){
            right = i;
            break;
        }
    }

    // find the left shoulder point
    for(int i = max_idx; i > 0; --i){
        if(Mn[i] >= shoulder && Mn[i-1] <= shoulder){
            left = i;
            break;
        }
    }

    if(left == -1 || right == -1){
//        delete [] Mn;
        free(Mn);
        if (max_idx*DIVIDE_FACTOR==0) {
            __android_log_print(ANDROID_LOG_VERBOSE, "debug", "out1 %d", max_idx * DIVIDE_FACTOR);
        }
        return max_idx*DIVIDE_FACTOR;
    }
    else{
        double middle = (right+left)*DIVIDE_FACTOR/2;
//        delete [] Mn;
        free(Mn);
        if (lround(middle)==0) {
            __android_log_print(ANDROID_LOG_VERBOSE, "debug", "out2 %d", lround(middle));
        }
        return lround(middle);
    }
}

// result, tx -> channel length = Ns
// rx -> recv_len = 8*(Ns+N0)
void channel_estimation_freq_multiple(double* result, double * tx, unsigned long int Ns, double * rx, unsigned long int recv_len, unsigned long int N0, int BW1, int BW2, int fs){
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "channel estimation %d %d %d",recv_len,Ns,N0);

    assert(  recv_len == 8*(Ns+N0) );
    short PN_seq[8] = {1, -1, -1, -1, -1, -1, 1, -1};

    fftw_complex *H_avg =(fftw_complex*)fftw_malloc(Ns*sizeof(fftw_complex));
    if(H_avg==NULL){fprintf(stderr,"malloc failed\n");exit(1);}

    for(int i = 0; i<8; ++i){
        fftw_complex *H_result =(fftw_complex*)fftw_malloc(Ns*sizeof(fftw_complex));
        if(H_result==NULL){fprintf(stderr,"malloc failed\n");exit(1);}

        channel_estimation_freq_single(H_result, tx, rx + i*(N0+Ns) + N0, Ns, BW1, BW2, fs);
        if(i == 0){
            for(unsigned long int j = 0; j < Ns; ++j ){
                H_avg[j][0] = PN_seq[i]*H_result[j][0];
                H_avg[j][1] = PN_seq[i]*H_result[j][1];
            }
        }
        else{
            for(unsigned long int j = 0; j < Ns; ++j ){
                H_avg[j][0] += PN_seq[i]*H_result[j][0];
                H_avg[j][1] += PN_seq[i]*H_result[j][1];
            }
        }

        fftw_free(H_result);
    }

    for(unsigned long int j = 0; j < Ns; ++j ){
        H_avg[j][0] /= (8*Ns); //H_result[j][0];
        H_avg[j][1] /= (8*Ns); //H_result[j][1];
    }

    fftw_complex *h=(fftw_complex*)fftw_malloc(Ns*sizeof(fftw_complex));
    if(h==NULL){fprintf(stderr,"malloc failed\n");exit(1);}

    fftw_plan ifft=fftw_plan_dft_1d(Ns, H_avg, h, FFTW_BACKWARD, FFTW_ESTIMATE);

    //double* abs_h = new double[Nu];

    fftw_execute(ifft);

    abs_complex(result, h, Ns);
//    char* out = getString_d(result,720);

    fftw_destroy_plan(ifft);

    fftw_free(H_avg);
    fftw_free(h);
}

JNIEXPORT void JNICALL
Java_com_example_nativeaudio_NativeAudio_naiserCorrTest(JNIEnv *env, jclass clazz, jdoubleArray tpre1,
                                                        jdoubleArray tpre2, jdoubleArray signal, jint corr_idx) {
    int bias=200;
    int win_size=2400;

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "naiser corr test");

    jint Nsignal1 = (*env)->GetArrayLength(env, signal);
    jint Ntx = (*env)->GetArrayLength(env, tpre1);
    jint Ntx2 = (*env)->GetArrayLength(env, tpre2);

    jdouble* signal1 = (*env)->GetDoubleArrayElements(env,signal, NULL);
    jdouble* pre1 = (*env)->GetDoubleArrayElements(env,tpre1, NULL);
    jdouble* pre2 = (*env)->GetDoubleArrayElements(env,tpre2, NULL);

    double* filtered=filter_d(signal1,Nsignal1);
//    char* filteredStr=getString_d(filtered,Nsignal1);

    // segment out after coarse xcorr
//    int corr_idx=182005;

    jint Nrx=Ntx2+(win_size*2);
    double* signal2 = calloc(Nrx,sizeof(double));
    memcpy(signal2,&filtered[corr_idx-win_size],Nrx*sizeof(double));

//    char* str1=getString_d(signal1,Nsignal1);
//    char* str2=getString_d(signal2,Nrx);

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "N %d",Nrx);

    int naiser_idx=naiser_corr(signal2, Nrx, 720, 480, 2, JNI_FALSE, .5, .8);
    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***naiser index %d",naiser_idx);

    double* signal3 = calloc(Ntx2,sizeof(double));
    memcpy(signal3,&signal2[naiser_idx+1-bias],Ntx2*sizeof(double));

    int BW1=1000;
    int BW2=5000;
    unsigned long int N0 = 480;
    double* h = calloc(Ntx,sizeof(double));

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "start h");

    channel_estimation_freq_multiple(h, pre1, Ntx, signal3, Ntx2, N0,  BW1,  BW2,  FS);

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "end h");

//    char* sig3Str = getString_d(h,Ntx);

    int max_idx = find_max(h, Ntx);
    int global_idx=(corr_idx-win_size+naiser_idx-bias)+max_idx;

    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***index %d %d %d %d",corr_idx,naiser_idx,max_idx,global_idx);
}

int findhpeak(double* h, int Ntx1, int bias) {
    double maxval=-1;
    int maxidx=-1;
    for (int i = 0; i < Ntx1; i++) {
        if (h[i] > maxval) {
            maxval=h[i];
            maxidx=i;
        }
    }
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug","h vals %f %f %f %f %f",h[0],h[1],h[2],h[3],h[4]);
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug","h vals %f %f %f %f %f",h[715],h[716],h[717],h[718],h[719]);
    __android_log_print(ANDROID_LOG_VERBOSE, "debug","h max %d %f",maxidx,maxval);

    int* peaksidxs = calloc(Ntx1,sizeof(int));
    int peakcounter=0;
    int prevpeakidx=0;
    int prevpeakval=0;
    for (int i = 1; i < maxidx+10; i++) {
        if (i < Ntx1 && h[i-1] < h[i] && h[i] > h[i+1] && h[i] > maxval*.65 && i-prevpeakidx >= 2) {
            peaksidxs[peakcounter++]=i;
        }
    }

    int midx_new=maxidx;
    for (int i = 0; i < peakcounter; i++) {
        if (maxidx-peaksidxs[i] <= bias) {
            midx_new=peaksidxs[i];
            break;
        }
    }

    free(peaksidxs);

    return midx_new;
}

// if output is positive, good it passes, else... fail
int corr2(int N, int xcorr_idx, double* filteredData, mycontext* cxt2, int globalOffset) {
    __android_log_print(ANDROID_LOG_VERBOSE, "debug2","corr2");
    int* out = calloc(2,sizeof(int));

    int outidx = -1;
    int offset=0;
    int Nu=cxt2->naiserTx1Count;
    int Ns2=cxt2->naiserTx2Count;
    unsigned long int N0 = cxt2->N0;
    jboolean CP = cxt2->CP;
    jint win_size = cxt2->win_size;
    jint bias = cxt2->bias;

//    jint Nrx=Ns2+(win_size*2);
    jint Nrx=Ns2+(win_size);
    int start_idx=xcorr_idx-win_size;
    int end_idx=start_idx+Nrx;

    __android_log_print(ANDROID_LOG_VERBOSE, "debug2","corr2 %d %d %d",start_idx,end_idx,N);
    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***corr2 cindex %d %d %d",outidx,Nu,Ns2);
//    if (xcorr_idx-win_size >= 0 && xcorr_idx+Ns2+win_size < N) {
    if (start_idx >= 0 && end_idx < N) {
        double* naiser_sig = calloc(Nrx,sizeof(double));

        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***naiser segment %d %d",start_idx,end_idx);
        memcpy(naiser_sig,&filteredData[start_idx],Nrx*sizeof(double));

//        char* str1=getString_d(naiser_sig,Nrx);
//        char* str2=getString_d(filteredData,N);

        clock_t begin = clock();

//        int* result = xcorr(filtered, refData2, N, N_ref2, 0,
//                            0, 4, 1.5, 960, .65,JNI_FALSE);
        int naiser_idx = naiser_corr(naiser_sig, Nrx, Nu, N0, 8, !CP, cxt2->naiserThresh, cxt2->naiserShoulder);

        clock_t end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        __android_log_print(ANDROID_LOG_VERBOSE, "time2", "%.4f", time_spent);

//        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***nindex %d",naiser_idx);
        if (naiser_idx>0) {
//            outidx = naiser_idx+(xcorr_idx-win_size);
//            outidx=xcorr_idx;
            outidx = xcorr_idx; // we don't care about the exact naiser_idx, only failure case
        }
        else {
            // did not pass naiser, false positive
//            free(naiser_sig);
            outidx = naiser_idx;
        }

        int start_idx2 = xcorr_idx - bias +1;
        int end_idx2 = start_idx2 + Ns2;
        if (start_idx2 >= 0 && end_idx2 < N && naiser_idx > 0) {
            int BW1=1000;
            int BW2=5000;

            int Ntx1=cxt2->naiserTx1Count;
            int Ntx2=cxt2->naiserTx2Count;
            double* h = calloc(Ntx1,sizeof(double));

            double* h_sig = calloc(Ns2,sizeof(double));
//            memcpy(h_sig,&naiser_sig[naiser_idx+1-bias],Ns2*sizeof(double));
//            memcpy(h_sig,&naiser_sig[win_size - bias +1],Ns2*sizeof(double));
            memcpy(h_sig,&filteredData[start_idx2],Ns2*sizeof(double));
//            char* str5=getString_d(h,Ntx1);

            __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***hindex %d %f %f %f %f %f",
                                xcorr_idx,filteredData[0],filteredData[1],filteredData[2],filteredData[3],filteredData[4]);
            __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***hindex %d %f %f %f %f %f",
                                xcorr_idx,h_sig[0],h_sig[1],h_sig[2],h_sig[3],h_sig[4]);

            begin = clock();

            channel_estimation_freq_multiple(h, cxt2->naiserTx1, Ntx1,
                                             h_sig, Ntx2,
                                             N0,  BW1,  BW2,  FS);
            int h_idx=findhpeak(h,Ntx1,cxt2->bias);

            end = clock();
            time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
            __android_log_print(ANDROID_LOG_VERBOSE, "time3", "%.4f", time_spent);

            free(h);
            free(h_sig);
            if (h_idx>0) {
                outidx = xcorr_idx - bias + h_idx;
            }
            __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***hindex %d %d",h_idx,outidx);
        }
        else {
            __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***h fail %d %d %d %d",naiser_idx, naiser_idx+1-bias, naiser_idx+Ns2-bias, N);
        }
        free(naiser_sig);
    }
    else {
        outidx = -1;
        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***nindex fail %d %d %d %d",xcorr_idx, xcorr_idx-win_size, xcorr_idx+Ns2+win_size, N);
    }

    return outidx;
}

int* xcorr_helper2(void* context, short* data, int globalOffset, int N) {
    mycontext* cxt = (mycontext*)context;
//    __android_log_print(ANDROID_LOG_VERBOSE, "debug2", "xcorr2_helper %d %d",globalOffset,N);

//    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***xcorr helper %d %d %d",
//                        naiser_index_to_process,cxt->processedSegments,cxt->getOneMoreFlag);

//    __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***normal corr");

    int N_ref = cxt->naiserTx2Count;

    double *refData = cxt->naiserTx2;

    double *filteredData = filter_s(data, N);
    for (int i = 0; i < N; i++) {
        filteredData[i] /= 10000.0;
    }
    double* filteredData2=filteredData+filt_offset;
    N -= filt_offset;

    clock_t begin = clock();
//        char* str1=getString_s(data,N);
//        char* str3=getString_d(filteredData,N);
//        char* str4=getString_d(refData,N_ref);

    int* xcorr_out = xcorr(filteredData2, refData, N, N_ref, cxt->queuedSegments,
                   globalOffset, cxt->xcorrthresh, cxt->minPeakDistance, cxt->seekback,
                   cxt->pthresh,cxt->getOneMoreFlag);
    //idx,val
    __android_log_print(ANDROID_LOG_VERBOSE, "debug6","before naiser %d %d %.0f", xcorr_out[0], xcorr_out[1],cxt->xcorrthresh);

    clock_t end = clock();

    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

//    __android_log_print(ANDROID_LOG_VERBOSE, "time1", "%.4f", time_spent);

    if (xcorr_out[0]<0) {
//        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***xcorr output neg %d", xcorr_out[0]);
    }
    else {
//        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***xcorr output pos %d", xcorr_out[0]);
    }

    int* result = calloc(3,sizeof(int));
    result[0] = xcorr_out[0]; //xcorr out
    result[1] = xcorr_out[1]; //xcorr val

    if (cxt->naiser && xcorr_out[0] > 0) {
        int idx = corr2(N,xcorr_out[0],filteredData2,cxt,globalOffset);
        result[2] = idx; // naiser out
//        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "***corr2 output %d %d %d",result[0],result[1],result[2]);
    }

    free(filteredData);

    return result;
}

JNIEXPORT void JNICALL
Java_com_example_nativeaudio_NativeAudio_testxcorr(JNIEnv *env, jclass clazz,jdoubleArray tempData,
                                                   jdoubleArray tempRefData1, jdoubleArray tempRefData2, jint N0,
                                                   jboolean CP) {
    jint N = (*env)->GetArrayLength(env, tempData);
    jint N_ref1 = (*env)->GetArrayLength(env, tempRefData1);
    jint N_ref2 = (*env)->GetArrayLength(env, tempRefData2);

    jdouble* data = (*env)->GetDoubleArrayElements(env,tempData, NULL);
    jdouble* refData1 = (*env)->GetDoubleArrayElements(env,tempRefData1, NULL);
    jdouble* refData2 = (*env)->GetDoubleArrayElements(env,tempRefData2, NULL);

    __android_log_print(ANDROID_LOG_VERBOSE, "debug","file lens %d %d %d",N,N_ref1,N_ref2);

    int filt_offset=127;
    for (int i = 0; i < 5; i++) {
        chirp_indexes[i]=-1;
    }
//    char* str1=getString_d(data,N);
//    char* str2=getString_d(refData2,N_ref2);

    double* filtered=filter_d(data,N);

    for (int i = 0; i < N; i++) {
        filtered[i]/=10000.0;
    }
    filtered=filtered+filt_offset;
    N -= filt_offset;

//    int* xcorr(double* filteredData, double* refData, int filtLength, int refLength, int i, int globalOffset, int xcorrthresh, double minPeakDistance, int seekback, double pthresh, jboolean getOneMoreFlag) {
    int* result = xcorr(filtered, refData2, N, N_ref2, 0,
                    0, 4, 1.5, 960, .65,JNI_FALSE);

//    char* str3=getString_d(filtered,N);

    int corr_idx= result[0];

    __android_log_print(ANDROID_LOG_VERBOSE, "debug","corrected xcorr after filter %d %d",result[0],result[1]);

    int win_size = 4800;
    int Ntx2= N_ref2;
    int Nu=N_ref1;
    int bias = 320;

    jint Nrx=Ntx2+(win_size*2);
    double* signal2 = calloc(Nrx,sizeof(double));

    int start_idx = corr_idx-win_size;
    int end_idx = start_idx+Nrx;
    if (start_idx < 0 || end_idx > N) {
        __android_log_print(ANDROID_LOG_VERBOSE, "debug","naiser bounds error %d %d %d",start_idx,end_idx,Nrx);
    }
    else {
        memcpy(signal2, &filtered[start_idx], Nrx * sizeof(double));

        char *nstr = getString_d(signal2, Nrx);
        int naiser_idx = naiser_corr(signal2, Nrx, Nu, N0, 8, !CP, .45, .8);
        __android_log_print(ANDROID_LOG_VERBOSE, "debug", "naiser %d", naiser_idx);

        if (naiser_idx >= 0) {
            double *h_sig = calloc(Ntx2, sizeof(double));
            int start_idx2 = corr_idx - bias + 1;
            int end_idx2 = start_idx2 + Ntx2;
            if (start_idx2 < 0 || end_idx2 > N) {
                __android_log_print(ANDROID_LOG_VERBOSE, "debug", "hsig bounds error %d %d %d",
                                    start_idx, end_idx, Nrx);
            } else {
                memcpy(h_sig, &filtered[start_idx2], Ntx2 * sizeof(double));

                //    char* out = getString_d(h_sig, 8640);

                double *h = calloc(Nu, sizeof(double));
                //    void channel_estimation_freq_multiple(double* result, double * tx, unsigned long int Ns, double * rx, unsigned long int recv_len, unsigned long int N0, int BW1, int BW2, int fs){

                channel_estimation_freq_multiple(h, refData1, Nu,
                                                 h_sig, Ntx2,
                                                 N0, 1000, 5000, 44100);
                //    char* out2 = getString_d(h, 720);

                int h_idx = findhpeak(h, Nu, bias);
                int global_idx = corr_idx - bias + h_idx;
                __android_log_print(ANDROID_LOG_VERBOSE, "debug", "h index %d", h_idx);
                __android_log_print(ANDROID_LOG_VERBOSE, "debug", "global index %d", global_idx);
            }
        }
    }
}

JNIEXPORT jint JNICALL
Java_com_example_nativeaudio_NativeAudio_getNumChirps(JNIEnv *env, jclass clazz) {
    return chirpsPlayed;
}