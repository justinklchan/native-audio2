package com.example.nativeaudio;

import android.content.Context;
import android.content.SharedPreferences;
import android.media.AudioManager;
import android.os.AsyncTask;
import android.os.CountDownTimer;
import android.preference.PreferenceManager;
import android.util.Log;
import android.widget.Button;
import android.widget.EditText;
import android.widget.Switch;
import android.widget.TextView;

import androidx.constraintlayout.widget.ConstraintLayout;
import androidx.core.widget.NestedScrollView;

import java.util.ArrayList;

public class Constants {
//    public static int bufferSize = 192*25; // 4800
//public static int bufferSize = 192*60; // 11520
//    public static int bufferSize = 192*250; // 48000
//    public static int bufferSize = 192*500; // 96000
    public static int minbuffersize;
    public static int bufferSize,bigBufferSize,bigBufferTimes;
    public static boolean stop=false;
    public static boolean initsleeping=false;
    static EditText et1,et2,et3,et4,et5,et6,et7,et8,et9,et10,et11,et12,et13,et14,et15;
    public static NestedScrollView sview;
    static int fs=44100;
    static float naiserThresh=0.1f, naiserShoulder=0.8f;
    static int win_size=4800;
    static int bias=320;
    static int initSleep=0;
    static float xcorrthresh=2f;
    public static float replyDelay = 3.0f;
//    static int calibSigLen = 4800;
    static boolean water,reply,naiser=true,runxcorr=true;
    static float vol=0.01f, minPeakDistance=1.5f;
    static int fileID=0;
    static int recTime=30;
    static int N0;
    static boolean CP;
    static TextView tv,tv2,distOut,debugPane,tv11;
    static Switch sw1,sw2,sw3,sw4;
    static Button recButton,stopButton,clearButton;
    static AsyncTask task;
    static CountDownTimer timer;
    static double[] pre1,pre2;
    static short[] sig=null;
    static int seekback=960;
    static float pthresh=.65f;
    static int rounds = 1;
    static float initialDelay = 3f;
    static float bufSizeInSeconds=.01f;
    static float bigBufferSizeInSeconds = .5f;
    static ConstraintLayout clayout;

    static ArrayList<Long> time_acc,time_gyro,time_acc_uncalib,time_gyro_uncalib,time_mag,time_mag_uncalib,time_pressure;
    static ArrayList<Float> accx,accy,accz;
    static ArrayList<Float> gyrox,gyroy,gyroz;
    static ArrayList<Float> accx_uncalib,accy_uncalib,accz_uncalib;
    static ArrayList<Float> gyrox_uncalib,gyroy_uncalib,gyroz_uncalib;
    static ArrayList<Float> magx,magy,magz;
    static ArrayList<Float> magx_uncalib,magy_uncalib,magz_uncalib;
    static ArrayList<Float> pressure_data;
    static boolean recordImu=false;
    static long tt;

    public static void init(Context cxt) {
        SharedPreferences prefs = PreferenceManager.getDefaultSharedPreferences(cxt);

        Constants.vol=prefs.getFloat("vol",vol);
        Constants.recTime=prefs.getInt("recTime",recTime);
        Constants.replyDelay=prefs.getFloat("replyDelay",replyDelay);
        Constants.initSleep=prefs.getInt("initSleep",initSleep);

        Constants.water=prefs.getBoolean("water",water);
        Constants.reply=prefs.getBoolean("reply",reply);
        Constants.naiser=prefs.getBoolean("naiser",naiser);
        Constants.xcorrthresh=prefs.getFloat("xcorrthresh",xcorrthresh);
        Constants.minPeakDistance=prefs.getFloat("minPeakDistance",minPeakDistance);
        Constants.fileID=prefs.getInt("fileID",fileID);
        Constants.naiserThresh=prefs.getFloat("naiserThresh",naiserThresh);
        Constants.naiserShoulder=prefs.getFloat("naiserShoulder",naiserShoulder);
        Constants.win_size=prefs.getInt("win_size",win_size);
        Constants.bias=prefs.getInt("bias",bias);
        Constants.seekback=prefs.getInt("seekback",seekback);
        Constants.pthresh=prefs.getFloat("pthresh",pthresh);
        Constants.rounds=prefs.getInt("rounds",rounds);
        Constants.runxcorr=prefs.getBoolean("runxcorr",runxcorr);
        Constants.initialDelay=prefs.getFloat("initialDelay",initialDelay);

        et1.setText(vol+"");
        et2.setText(recTime+"");
        et3.setText(replyDelay+"");
        et4.setText(initSleep+"");
        et5.setText(xcorrthresh+"");
        et6.setText(minPeakDistance+"");
        et7.setText(fileID+"");
        et8.setText(naiserThresh+"");
        et9.setText(naiserShoulder+"");
        et10.setText(win_size+"");
        et11.setText(bias+"");
        et12.setText(seekback+"");
        et13.setText(pthresh+"");
        et14.setText(rounds+"");
        et15.setText(initialDelay+"");
        sw1.setChecked(water);
        sw2.setChecked(reply);
        sw3.setChecked(naiser);
        sw4.setChecked(runxcorr);

        AudioManager am = (AudioManager) cxt.getSystemService(Context.AUDIO_SERVICE);
        minbuffersize=Integer.parseInt(am.getProperty(AudioManager.PROPERTY_OUTPUT_FRAMES_PER_BUFFER));

        if (sw1.isChecked()) {
            Constants.sw1.setText("Water");
        }
        else {
            Constants.sw1.setText("Air");
        }

        float bufSizeInSamples = fs*bufSizeInSeconds;
        bufferSize=((int)Math.ceil(bufSizeInSamples/minbuffersize))*minbuffersize;
//        bufferSize+=(minbuffersize*3);

        float bigBufferSizeInSamples = fs*bigBufferSizeInSeconds;
        bigBufferTimes=((int)Math.ceil(bigBufferSizeInSamples/bufferSize));
        bigBufferSize = bigBufferTimes*bufferSize;
        Log.e("asdf","BUFFER "+minbuffersize+","+bufferSize+","+bigBufferSize);
        loadData(cxt);
    }

    public static void toggleUI() {
        Constants.recButton.setEnabled(!Constants.recButton.isEnabled());
        Constants.stopButton.setEnabled(!Constants.stopButton.isEnabled());
        Constants.clearButton.setEnabled(!Constants.clearButton.isEnabled());
        Constants.sw1.setEnabled(!Constants.sw1.isEnabled());
        Constants.sw2.setEnabled(!Constants.sw2.isEnabled());
        Constants.sw3.setEnabled(!Constants.sw3.isEnabled());
        Constants.sw4.setEnabled(!Constants.sw4.isEnabled());
        Constants.et1.setEnabled(!Constants.et1.isEnabled());
        Constants.et2.setEnabled(!Constants.et2.isEnabled());
        Constants.et3.setEnabled(!Constants.et3.isEnabled());
        Constants.et4.setEnabled(!Constants.et4.isEnabled());
        Constants.et5.setEnabled(!Constants.et5.isEnabled());
        Constants.et6.setEnabled(!Constants.et6.isEnabled());
        Constants.et7.setEnabled(!Constants.et7.isEnabled());
        Constants.et8.setEnabled(!Constants.et8.isEnabled());
        Constants.et9.setEnabled(!Constants.et9.isEnabled());
        Constants.et10.setEnabled(!Constants.et10.isEnabled());
        Constants.et11.setEnabled(!Constants.et11.isEnabled());
        Constants.et12.setEnabled(!Constants.et12.isEnabled());
        Constants.et13.setEnabled(!Constants.et13.isEnabled());
        Constants.et14.setEnabled(!Constants.et14.isEnabled());
        Constants.et15.setEnabled(!Constants.et15.isEnabled());
    }

    public static void loadData(Context cxt) {
        String txt = "null";
        if (Constants.fileID==0) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.online1);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.train_sig1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.train_sig2));
            txt="online1";
            N0=480;
        }
        else if (Constants.fileID==1) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_720_360_cp);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.cp_train_sig1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.cp_train_sig2));
            txt="signal_N_720_360_CP";
            N0=360;
            CP=true;
        }
        else if (Constants.fileID==2) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_720_360_gi);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.gi_train_sig1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.gi_train_sig2));
            txt="signal_N_720_360_GI";
            N0=360;
            CP=false;
        }
        else if (Constants.fileID==3) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_720_360_half);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.half_train_sig1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.half_train_sig2));
            txt="signal_N_720_360_HALF";
            N0=360;
            CP=false;
        }
        else if (Constants.fileID==4) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.n_1080_360_half_signal);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n_1080_360_half_signal_train_sig1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n_1080_360_half_signal_train_sig2));
            txt="signal_N_1080_360_HALF";
            N0=360;
            CP=false;
        }
        else if (Constants.fileID==5) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.n_1260_480_half_signal);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n_1260_480_half_train_sig1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n_1260_480_half_train_sig2));
            txt="signal_N_1260_480_HALF";
            N0=480;
            CP=false;
        }
        else if (Constants.fileID==6) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_960_360_1000_3000);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n960_360_1000_3000_t1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n960_360_1000_3000_t2));
            txt="n960_360_1000_3000_t1";
            N0=360;
            CP=false;
        }
        else if (Constants.fileID==7) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_960_360_1000_5000);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n960_360_1000_5000_t1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n960_360_1000_5000_t2));
            txt="n960_360_1000_5000";
            N0=360;
            CP=false;
        }
        else if (Constants.fileID==8) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_960_360_1000_9000);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n960_360_1000_9000_t1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n960_360_1000_9000_t2));
            txt="n960_360_1000_9000";
            N0=360;
            CP=false;
        }
        else if (Constants.fileID==9) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_1260_480_1000_3000);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n1260_480_1000_3000_t1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n1260_480_1000_3000_t2));
            txt="n1260_480_1000_3000";
            N0=480;
            CP=false;
        }
        else if (Constants.fileID==10) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_1260_480_1000_5000);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n1260_480_1000_5000_t1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n1260_480_1000_5000_t2));
            txt="n1260_480_1000_5000";
            N0=480;
            CP=false;
        }
        else if (Constants.fileID==11) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_1260_480_1000_6000);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n1260_480_1000_6000_t1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n1260_480_1000_6000_t2));
            txt="n1260_480_1000_6000";
            N0=480;
            CP=false;
        }
        else if (Constants.fileID==12) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_1260_480_1000_7000);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n1260_480_1000_7000_t1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n1260_480_1000_7000_t2));
            txt="n1260_480_1000_7000";
            N0=480;
            CP=false;
        }
        else if (Constants.fileID==13) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_1260_480_1000_9000);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n1260_480_1000_9000_t1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n1260_480_1000_9000_t2));
            txt="n1260_480_1000_9000";
            N0=480;
            CP=false;
        }
        else if (Constants.fileID==14) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_2160_480_1000_3000);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n2160_480_1000_3000_t1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n2160_480_1000_3000_t2));
            txt="n1260_480_1000_3000";
            N0=480;
            CP=false;
        }
        else if (Constants.fileID==15) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_2160_480_1000_5000);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n2160_480_1000_5000_t1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n2160_480_1000_5000_t2));
            txt="n1260_480_1000_5000";
            N0=480;
            CP=false;
        }
        else if (Constants.fileID==16) {
            sig=FileOperations.readrawasset_binary(cxt, R.raw.signal_2160_480_1000_9000);
            pre1=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n2160_480_1000_9000_t1));
            pre2=Utils.convert(FileOperations.readrawasset_binary(cxt,R.raw.n2160_480_1000_9000_t2));
            txt="n1260_480_1000_9000";
            N0=480;
            CP=false;
        }

        if (Utils.max(pre1) > 2) {
            Utils.div(pre1, 32767);
        }
        if (Utils.max(pre2) > 2) {
            Utils.div(pre2,32767);
        }

        String finalTxt = txt;
        (NativeAudio.av).runOnUiThread(new Runnable() {
            @Override
            public void run() {
                tv11.setText(finalTxt);
            }
        });
        Log.e("asdf","load file "+sig.length);
    }
}
