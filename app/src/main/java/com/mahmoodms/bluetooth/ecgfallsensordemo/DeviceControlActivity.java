package com.mahmoodms.bluetooth.ecgfallsensordemo;

import android.app.ActionBar;
import android.app.Activity;
import android.bluetooth.BluetoothDevice;
import android.bluetooth.BluetoothGatt;
import android.bluetooth.BluetoothGattCharacteristic;
import android.bluetooth.BluetoothGattDescriptor;
import android.bluetooth.BluetoothGattService;
import android.bluetooth.BluetoothManager;
import android.bluetooth.BluetoothProfile;
import android.content.BroadcastReceiver;
import android.content.Context;
import android.content.Intent;
import android.content.IntentFilter;
import android.content.pm.ActivityInfo;
import android.graphics.Color;
import android.graphics.Typeface;
import android.graphics.drawable.ColorDrawable;
import android.net.Uri;
import android.os.BatteryManager;
import android.os.Build;
import android.os.Bundle;
import android.os.Environment;
import android.os.Handler;
import android.provider.Settings;
import android.support.v4.app.NavUtils;
import android.telephony.SmsManager;
import android.util.Log;
import android.view.Menu;
import android.view.MenuItem;
import android.view.View;
import android.view.WindowManager;
import android.widget.Button;
import android.widget.CompoundButton;
import android.widget.Switch;
import android.widget.TextView;
import android.widget.Toast;

import com.androidplot.Plot;
import com.androidplot.util.Redrawer;
import com.androidplot.xy.BoundaryMode;
import com.androidplot.xy.LineAndPointFormatter;
import com.androidplot.xy.SimpleXYSeries;
import com.androidplot.xy.XYPlot;
import com.androidplot.xy.XYStepMode;
import com.beele.BluetoothLe;
import com.opencsv.CSVWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 * Created by mahmoodms on 5/31/2016.
 */

public class DeviceControlActivity extends Activity implements BluetoothLe.BluetoothLeListener {
    private final static String TAG = DeviceControlActivity.class.getSimpleName();
    //LocalVars
    private String mDeviceName;
    private String mDeviceAddress;
    private String mDeviceAddress2;
    private boolean mConnected;
    private boolean mConnected2;
    //Class instance variable
    private BluetoothLe mBluetoothLe;
    private BluetoothManager mBluetoothManager = null;
    private BluetoothGatt mBluetoothGatt = null;
    private BluetoothGatt mBluetoothGatt2 = null;
    private BluetoothDevice mBluetoothDevice;
    private BluetoothDevice mBluetoothDevice2;
    //Layout - TextViews and Buttons

    private TextView mEegValsTextView;
    private TextView mBatteryLevel;
    private TextView mDataRate;
    private Button mExportButton;
    private Switch mFilterSwitch;

    private long mLastTime;
    private long mLastTime2;
    private long mCurrentTime;
    private long mCurrentTime2;
    private int mTimeSeconds = 0;
    private String mLastRssi;

    private ScheduledThreadPoolExecutor executor = new ScheduledThreadPoolExecutor(1);

    private boolean filterData = false;
    private int dataCount = 0;
    private int points = 0;
    private Menu menu;
    //RSSI Stuff:
    private static final int RSSI_UPDATE_TIME_INTERVAL = 2000;
    private Handler mTimerHandler = new Handler();

    private boolean mTimerEnabled = false;
    /**
     * Initialize Plot:
     */
    private XYPlot eegPlot;
    private Redrawer redrawer;
    private SimpleXYSeries eegDataSeries1;
    private static final int HISTORY_SIZE = 1000;
    private static final int HISTORY_SECONDS = 4;
    private boolean plotImplicitXVals = false;
    private int DATA_RATE_SAMPLES_PER_SECOND = 0;
    //for bounds:
    private BoundaryMode currentBM = BoundaryMode.AUTO;
    //Data Variables:
    private int batteryWarning = 20;//%
    final private boolean terminate = true;
    final private boolean initialize = false;
    private String fileTimeStamp = "";
    private long periodShort = 3500; // 3.5ms (changed from 3.7ms on 10/11 @  7:50pm)
    private int eegIndex = 0;
    private double dataRate;
    private double dataRate2;
    /* Notification stuff */
    private Button mC1Button;
    private Button mC0Button;
    private Button mC2Button;
    private Button mC3Button;
    private Button mC4Button;
    private BluetoothGattService mLedService = null;

    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_device_control);
        //Set orientation of device based on screen type/size:
        setRequestedOrientation(ActivityInfo.SCREEN_ORIENTATION_LANDSCAPE);
        //Recieve Intents:
        Intent intent = getIntent();
        mDeviceName = intent.getStringExtra(AppConstant.EXTRAS_DEVICE_NAME);
        mDeviceAddress = intent.getStringExtra(AppConstant.EXTRAS_DEVICE_ADDRESS);
        mDeviceAddress2 = "C7:93:2D:DD:64:D1";
        Log.d(TAG, "mDeviceAddress(1): "+mDeviceAddress);
        Log.d(TAG, "mDeviceAddress(2): "+mDeviceAddress2);
        //Set up action bar:
        if (getActionBar() != null) {
            getActionBar().setDisplayHomeAsUpEnabled(true);
        }
        ActionBar actionBar = getActionBar();
        actionBar.setBackgroundDrawable(new ColorDrawable(Color.parseColor("#6078ef")));

        //Flag to keep screen on (stay-awake):
        getWindow().addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);
        //Set up TextViews
        mEegValsTextView = (TextView) findViewById(R.id.ecgValue);
        mExportButton = (Button) findViewById(R.id.button_export);
        mFilterSwitch = (Switch) findViewById(R.id.filterSwitch);
        mBatteryLevel = (TextView) findViewById(R.id.batteryText);
        mDataRate = (TextView) findViewById(R.id.dataRate);
        mDataRate.setText("...");
        //Initialize Bluetooth
        ActionBar ab = getActionBar();
        ab.setTitle(mDeviceName);
        ab.setSubtitle(mDeviceAddress);
        initializeBluetooth();
        // Initialize our XYPlot reference:
        eegDataSeries1 = new SimpleXYSeries("EEG Data Ch 1 (V)");
        eegPlot = (XYPlot) findViewById(R.id.eegPlot);
        //Todo: Graph temporarily uses data index - find alternative to implicit XVals→(seconds)
        if (plotImplicitXVals) {
            eegDataSeries1.useImplicitXVals();
            eegPlot.setDomainBoundaries(0, HISTORY_SIZE, BoundaryMode.FIXED);
            eegPlot.setDomainStepMode(XYStepMode.INCREMENT_BY_VAL);
            eegPlot.setDomainStepValue(HISTORY_SIZE/5);
        } else {
            eegPlot.setDomainBoundaries(0, HISTORY_SECONDS, BoundaryMode.FIXED);
            eegPlot.setDomainStepMode(XYStepMode.INCREMENT_BY_VAL);
            eegPlot.setDomainStepValue(HISTORY_SECONDS / 4);
        }
        if(filterData) {
            eegPlot.setRangeBoundaries(-2.5, 2.5, BoundaryMode.AUTO); //EMG only!
            eegPlot.setRangeStepValue(1);
        } /*else {
            eegPlot.setRangeBoundaries(-2.5, 2.5, BoundaryMode.FIXED);
            eegPlot.setRangeStepValue(0.5);
        }*/
        eegPlot.setRangeStepMode(XYStepMode.INCREMENT_BY_VAL);
        eegPlot.setDomainLabel("Time (seconds)");
        eegPlot.getDomainLabelWidget().pack();
        eegPlot.setRangeLabel("Voltage (mV)");
        eegPlot.getRangeLabelWidget().pack();
        eegPlot.setRangeValueFormat(new DecimalFormat("#.###"));
        eegPlot.setDomainValueFormat(new DecimalFormat("#"));
        eegPlot.getDomainLabelWidget().getLabelPaint().setColor(Color.BLACK);
        eegPlot.getDomainLabelWidget().getLabelPaint().setTextSize(20);
        eegPlot.getRangeLabelWidget().getLabelPaint().setColor(Color.BLACK);
        eegPlot.getRangeLabelWidget().getLabelPaint().setTextSize(20);
        eegPlot.getGraphWidget().getDomainTickLabelPaint().setColor(Color.BLACK);
        eegPlot.getGraphWidget().getRangeTickLabelPaint().setColor(Color.BLACK);
        eegPlot.getGraphWidget().getDomainTickLabelPaint().setTextSize(23); //was 36
        eegPlot.getGraphWidget().getRangeTickLabelPaint().setTextSize(23);
        eegPlot.getGraphWidget().getDomainGridLinePaint().setColor(Color.WHITE);
        eegPlot.getGraphWidget().getRangeGridLinePaint().setColor(Color.WHITE);
        eegPlot.getLegendWidget().getTextPaint().setColor(Color.BLACK);
        eegPlot.getLegendWidget().getTextPaint().setTextSize(20);
        eegPlot.getTitleWidget().getLabelPaint().setTextSize(20);
        eegPlot.getTitleWidget().getLabelPaint().setColor(Color.BLACK);
        LineAndPointFormatter lineAndPointFormatter1 = new LineAndPointFormatter(Color.BLACK, null, null, null);
        lineAndPointFormatter1.getLinePaint().setStrokeWidth(3);
        LineAndPointFormatter lineAndPointFormatter2 = new LineAndPointFormatter(Color.RED, null, null, null);
        lineAndPointFormatter2.getLinePaint().setStrokeWidth(4);
        LineAndPointFormatter lineAndPointFormatter3 = new LineAndPointFormatter(Color.BLUE, null, null, null);
        lineAndPointFormatter3.getLinePaint().setStrokeWidth(2);
        LineAndPointFormatter lineAndPointFormatter4 = new LineAndPointFormatter(Color.parseColor("#19B52C"), null, null, null);
        lineAndPointFormatter4.getLinePaint().setStrokeWidth(3);
        eegPlot.addSeries(eegDataSeries1, lineAndPointFormatter1);
        redrawer = new Redrawer(
                Arrays.asList(new Plot[]{eegPlot}),
                100, false);
        mExportButton.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                try {
                    exportFile(false, true, "", 0.0);
                } catch (IOException e) {
                    Log.e("IOException", e.toString());
                }
            }
        });
        eegIndex = 750;
        makeFilterSwitchVisible(false);
        mFilterSwitch.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                filterData = isChecked;
                if(filterData) {
                    if(eegDataSeries1.size()>0) {
                        double average = 0;
                        for (int i = 0; i < eegDataSeries1.size(); i++) {
                            average+=(double)eegDataSeries1.getY(i);
                            average/=(eegDataSeries1.size());
                        }
                        eegPlot.setRangeBoundaries(average-0.010, average+0.010, BoundaryMode.AUTO);
                        eegPlot.setRangeStepValue(0.4);
                    }
                } else {
                    if(eegDataSeries1.size()>0) {
                        double average = 0;
                        if(eegDataSeries1.size()>250) {
                            for (int i = 0; i < 250; i++) {
                                average+=(double)eegDataSeries1.getY(i);
                                average/=(250.0);
                            }
                        } else {
                            for (int i = 0; i < eegDataSeries1.size(); i++) {
                                average+=(double)eegDataSeries1.getY(i);
                                average/=(eegDataSeries1.size());
                            }
                        }
                        eegPlot.setRangeBoundaries(average-0.010, average+0.010, BoundaryMode.AUTO);
                        eegPlot.setRangeStepValue(0.1);
                    }
                }
            }
        });
        mLastTime = System.currentTimeMillis();
        mLastTime2 = mLastTime;
        Arrays.fill(ECGBufferUnfiltered, 0.0);
        this.registerReceiver(this.mBatInfoReceiver, new IntentFilter(Intent.ACTION_BATTERY_CHANGED));
        if(mDeviceName.equals("ECG 1000Hz")) {
            DATA_RATE_SAMPLES_PER_SECOND = 1000;
        } else if (mDeviceName.equals("ECG 500Hz")||mDeviceName.equals("EEG 500Hz")) {
            DATA_RATE_SAMPLES_PER_SECOND = 500;
        } else if (mDeviceName.equals("ECG 250Hz")||mDeviceName.equals("ECGSensor")||mDeviceName.equals("EEG 250Hz")) {
            DATA_RATE_SAMPLES_PER_SECOND = 250;
        } else {
            DATA_RATE_SAMPLES_PER_SECOND = 250;
        }
        //TODO: Test fHC with dummy data:
        double[] ch1 = new double[1000];
        double[] ch2 = new double[1000];
        double[] ch3 = new double[1000];
        double[] ch4 = new double[1000];
        boolean eogOnly = false;
        Arrays.fill(ch1,0.0);
        Arrays.fill(ch2,0.0);
        Arrays.fill(ch3,0.0);
        Arrays.fill(ch4,0.0);
        double[] yfit = jfullHybridClassifier(ch1, ch2, ch3, ch4, eogOnly);
        Log.e(TAG, "fHC TEST ARRAY: "+Arrays.toString(yfit));
        mC0Button = (Button) findViewById(R.id.c0_off);
        mC1Button = (Button) findViewById(R.id.c1button);
        mC2Button = (Button) findViewById(R.id.c2button);
        mC3Button = (Button) findViewById(R.id.c3button);
        mC4Button = (Button) findViewById(R.id.c4button);
        mC0Button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                if(mConnected) {
                    byte[] bytes = new byte[1];
                    bytes[0] = (byte) 0x00;
                    if(mLedService!=null)
                    mBluetoothLe.writeCharacteristic(mBluetoothGatt2, mLedService.getCharacteristic(AppConstant.CHAR_WHEELCHAIR_CONTROL),bytes);
                }
            }
        });
        mC1Button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                if(mConnected) {
                    byte[] bytes = new byte[1];
                    bytes[0] = (byte) 0x01;
                    if(mLedService!=null)
                    mBluetoothLe.writeCharacteristic(mBluetoothGatt2, mLedService.getCharacteristic(AppConstant.CHAR_WHEELCHAIR_CONTROL),bytes);
                }
            }
        });
        mC2Button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                if(mConnected) {
                    byte[] bytes = new byte[1];
                    bytes[0] = (byte) 0x0F;
                    if(mLedService!=null)
                    mBluetoothLe.writeCharacteristic(mBluetoothGatt2, mLedService.getCharacteristic(AppConstant.CHAR_WHEELCHAIR_CONTROL),bytes);
                }
            }
        });
        mC3Button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                if(mConnected) {
                    byte[] bytes = new byte[1];
                    bytes[0] = (byte) 0xF0;
                    if(mLedService!=null)
                    mBluetoothLe.writeCharacteristic(mBluetoothGatt2, mLedService.getCharacteristic(AppConstant.CHAR_WHEELCHAIR_CONTROL),bytes);
                }
            }
        });
        mC4Button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                if(mConnected) {
                    byte[] bytes = new byte[1];
                    bytes[0] = (byte) 0xFF;
                    if(mLedService!=null)
                    mBluetoothLe.writeCharacteristic(mBluetoothGatt2, mLedService.getCharacteristic(AppConstant.CHAR_WHEELCHAIR_CONTROL),bytes);
                }
            }
        });
    }

    private String androidDeviceBatteryLevel = "-1%";
    private String androidDeviceBatteryStatus = "-1";
    private BroadcastReceiver mBatInfoReceiver = new BroadcastReceiver(){
        @Override
        public void onReceive(Context ctxt, Intent intent) {
            int level = intent.getIntExtra(BatteryManager.EXTRA_LEVEL, 0);
            int status = intent.getIntExtra(BatteryManager.EXTRA_STATUS, -1);
            androidDeviceBatteryLevel = String.valueOf(level) + "%";
            if(status == BatteryManager.BATTERY_STATUS_CHARGING || status == BatteryManager.BATTERY_STATUS_FULL) {
                androidDeviceBatteryStatus = "Charging";
            } else {
                androidDeviceBatteryStatus = "Discharging";
            }
        }
    };

    public String getTimeStamp() {
        return new SimpleDateFormat("yyyy.MM.dd_HH.mm.ss").format(new Date());
    }
    public String getTimeStamp2() {
        return new SimpleDateFormat("MM/dd/yyyy - HH:mm:ss").format(new Date())+" EST\r\n";
    }

    //Write Log File:
    private boolean fileLogInitialized = false;
    private File logFile;
    public void exportLogFile(boolean init, String dataToWrite) {
        if(init) {
            Log.i("exportLogFile", "generated Log file");
            root = Environment.getExternalStorageDirectory();
            File dir = new File(root+"/EEGDataLogs");
            boolean mkdirsA = dir.mkdirs();
            logFile = new File(dir,"Log_"+getTimeStamp()+".txt");
            if (!logFile.exists()) {
                try {
                    logFile.createNewFile();
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }
            fileLogInitialized = true;
        } else {
            if(logFile.exists()) {
                try {
                    BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(logFile,true));
                    bufferedWriter.append(dataToWrite);
                    bufferedWriter.newLine();
                    bufferedWriter.close();
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }
        }
    }

    private boolean fileExportInitialized = false;
    private CSVWriter csvWriter;
    private File file;
    private File root;
    private String[] valueCsvWrite = new String[1];
    private String[] valueCsvWrite2 = new String[4];
    private int exportFileDataPointCounter = 0;
    private int exportFilePart = 1;
    public void exportFile(boolean init, boolean terminateExport,
                           String fileName, double eegData1) throws IOException {
        if (init) {
            root = Environment.getExternalStorageDirectory();
            fileTimeStamp = fileName;
            fileExportInitialized = true;
        } else {
            if (fileTimeStamp == null || fileTimeStamp.equals("") || !fileExportInitialized) {
                fileTimeStamp = "EEGSensorData_" + getTimeStamp();
            }
        }
        if (root.canWrite() && init) {
            File dir = new File(root.getAbsolutePath() + "/DataDirectory");
            boolean mkdirsA = dir.mkdirs();
            file = new File(dir, fileTimeStamp + "_part"+ String.valueOf(exportFilePart) + ".csv")
            ;
            csvWriter = new CSVWriter(new FileWriter(file));
            Log.d("New File Generated", fileTimeStamp + "_part"+ String.valueOf(exportFilePart) + ".csv");
            exportLogFile(false, "NEW FILE GENERATED: "+fileTimeStamp + "_part"+ String.valueOf(exportFilePart) + ".csv\r\n\r\n");
            if(exportFilePart!=1)exportLogFile(false, getDetails());
        }
        //Write Data to File (if init & terminateExport are both false)
        if (!init && !terminateExport) {
            if(exportFileDataPointCounter<1048575) {
                valueCsvWrite[0] = eegData1 + "," + eegData1 + "," + eegData1 + "," + eegData1 + "";
                csvWriter.writeNext(valueCsvWrite);
                exportFileDataPointCounter++;
            } else {
                valueCsvWrite[0] = eegData1 + "";
                csvWriter.writeNext(valueCsvWrite);
                csvWriter.flush();
                csvWriter.close();
                exportFileDataPointCounter=0;
                exportFilePart++;
                //generate new file:
                exportFile(true, false, fileTimeStamp,0);
            }
        }
        if (terminateExport) {
            csvWriter.flush();
            csvWriter.close();
            Uri uii;
            uii = Uri.fromFile(file);
            Intent exportData = new Intent(Intent.ACTION_SEND);
            exportData.putExtra(Intent.EXTRA_SUBJECT, "ECG Data Export Details");
            exportData.putExtra(Intent.EXTRA_STREAM, uii);
            exportData.setType("text/html");
            startActivity(exportData);
        }

    }

    public void exportFile(boolean init, boolean terminateExport,
                           String fileName, double eegData1, double eegData2, double eegData3, double eegData4) throws IOException {
        if (!init && !terminateExport) {
            if(exportFileDataPointCounter<1048575) {
                valueCsvWrite2[0] = eegData1 + "";
                valueCsvWrite2[1] = eegData2 + "";
                valueCsvWrite2[2] = eegData3 + "";
                valueCsvWrite2[3] = eegData4 + "";
                csvWriter.writeNext(valueCsvWrite2,false);
                exportFileDataPointCounter++;
            } else {
                valueCsvWrite[0] = eegData1 + "";
                csvWriter.writeNext(valueCsvWrite);
                csvWriter.flush();
                csvWriter.close();
                exportFileDataPointCounter=0;
                exportFilePart++;
                //generate new file:
                exportFile(true, false, fileTimeStamp,0);
            }
        }
    }

    @Override
    public void onResume() {
        makeFilterSwitchVisible(true);
        int fHCMain = jmainFHC(initialize);
        Log.i("fHCMain","INITIALIZED: "+String.valueOf(fHCMain));
//        int feegcfilt = jmainEegFilt(initialize);
        double[] array0 = new double[1000];
        Arrays.fill(array0,0.0);
        double[] retArray = jeegcfilt(array0);
        String fileTimeStampConcat = "EEGSensorData_" + getTimeStamp();
        Log.d("onResume-timeStamp", fileTimeStampConcat);
        //TODO (IF ECG/EMG PRESENT ONLY!!!) → We're creating a lot of empty files!!!
        if(!fileLogInitialized) {
            exportLogFile(true, "");
        }
        if(!fileExportInitialized) {
            try {
                exportFile(true, false, fileTimeStampConcat, 0.0);
            } catch (IOException ex) {
                Log.e("IOEXCEPTION:", ex.toString());
            }
        }
        redrawer.start();
        super.onResume();
    }
    
    @Override
    protected void onPause() {
        redrawer.pause();
        makeFilterSwitchVisible(false);
        super.onPause();
    }

    private void initializeBluetooth() {
        mBluetoothManager = (BluetoothManager) getSystemService(Context.BLUETOOTH_SERVICE);
        mBluetoothDevice = mBluetoothManager.getAdapter().getRemoteDevice(mDeviceAddress);
        //TODO: Connect To Another Device.
        mBluetoothDevice2 = mBluetoothManager.getAdapter().getRemoteDevice(mDeviceAddress2);
        mBluetoothLe = new BluetoothLe(this, mBluetoothManager, this);
        mBluetoothGatt = mBluetoothLe.connect(mBluetoothDevice, false);
        mBluetoothGatt2 = mBluetoothLe.connect(mBluetoothDevice2, false);
    }

    private void setNameAddress(String name_action, String address_action) {
        MenuItem name = menu.findItem(R.id.action_title);
        MenuItem address = menu.findItem(R.id.action_address);
        name.setTitle(name_action);
        address.setTitle(address_action);
        invalidateOptionsMenu();
    }

    @Override
    protected void onDestroy() {
        redrawer.finish();
        mBluetoothLe.disconnect(mBluetoothGatt);
        mBluetoothLe.disconnect(mBluetoothGatt2);
        this.unregisterReceiver(mBatInfoReceiver);
        stopMonitoringRssiValue();
        super.onDestroy();
    }

    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
        getMenuInflater().inflate(R.menu.menu_device_control, menu);
        getMenuInflater().inflate(R.menu.actionbar_item, menu);
        if (mConnected) {
            menu.findItem(R.id.menu_connect).setVisible(false);
            menu.findItem(R.id.menu_disconnect).setVisible(true);
        } else {
            menu.findItem(R.id.menu_connect).setVisible(true);
            menu.findItem(R.id.menu_disconnect).setVisible(false);
        }
        this.menu = menu;
        setNameAddress(mDeviceName, mDeviceAddress);
        return true;
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        switch (item.getItemId()) {
            case R.id.menu_connect:
                if (mBluetoothLe != null)
                    mBluetoothLe.connect(mBluetoothDevice, false);
//                    mBluetoothLe.connect(mBluetoothDevice2, false);
                connect();
                return true;
            case R.id.menu_disconnect:
                if (mBluetoothLe != null) {
                    if (mBluetoothGatt != null) {
                        mBluetoothLe.disconnect(mBluetoothGatt);
                    }
                    if(mBluetoothGatt2 != null) {
                        mBluetoothLe.disconnect(mBluetoothGatt2);
                    }
                }
                return true;
            case android.R.id.home:
                if (mBluetoothLe != null) {
                    if(mBluetoothGatt!=null) {
                        mBluetoothLe.disconnect(mBluetoothGatt);
                    }
                    if(mBluetoothGatt2 != null) {
                        mBluetoothLe.disconnect(mBluetoothGatt2);
                    }
                }
                NavUtils.navigateUpFromSameTask(this);
                onBackPressed();
                return true;
        }
        return super.onOptionsItemSelected(item);
    }

    private void connect() {
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                MenuItem menuItem = menu.findItem(R.id.action_status);
                menuItem.setTitle("Connecting...");
            }
        });
    }

    @Override
    public void onServicesDiscovered(BluetoothGatt gatt, int status) {
        Log.i(TAG, "onServicesDiscovered");
        if (status == BluetoothGatt.GATT_SUCCESS) {
            for (BluetoothGattService service : gatt.getServices()) {
                if ((service == null) || (service.getUuid() == null)) {
                    continue;
                }
                if (AppConstant.SERVICE_DEVICE_INFO.equals(service.getUuid())) {
                    //Read the device serial number
                    mBluetoothLe.readCharacteristic(gatt, service.getCharacteristic(AppConstant.CHAR_SERIAL_NUMBER));
                    //Read the device software version
                    mBluetoothLe.readCharacteristic(gatt, service.getCharacteristic(AppConstant.CHAR_SOFTWARE_REV));
                }
                if(AppConstant.SERVICE_WHEELCHAIR_CONTROL.equals(service.getUuid())) {
                    mLedService = service;
                    Log.i(TAG,"BLE Wheelchair Control Service found");
                }

                if(AppConstant.SERVICE_EEG_SIGNAL.equals(service.getUuid())) {
                    makeFilterSwitchVisible(true);
                    mBluetoothLe.setCharacteristicNotification(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CH1_SIGNAL), true);
                    mBluetoothLe.setCharacteristicNotification(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CH2_SIGNAL), true);
                    mBluetoothLe.setCharacteristicNotification(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CH3_SIGNAL), true);
                    mBluetoothLe.setCharacteristicNotification(gatt, service.getCharacteristic(AppConstant.CHAR_EEG_CH4_SIGNAL), true);
                }

                if (AppConstant.SERVICE_BATTERY_LEVEL.equals(service.getUuid())) {
                    //Read the device battery percentage
                    mBluetoothLe.readCharacteristic(gatt, service.getCharacteristic(AppConstant.CHAR_BATTERY_LEVEL));
                    mBluetoothLe.setCharacteristicNotification(gatt, service.getCharacteristic(AppConstant.CHAR_BATTERY_LEVEL), true);
                }
                if (AppConstant.SERVICE_MPU.equals(service.getUuid())) {
                    mBluetoothLe.setCharacteristicNotification(gatt, service.getCharacteristic(AppConstant.CHAR_MPU_COMBINED), true);
                }
            }
        }
    }

    private void makeFilterSwitchVisible(final boolean visible) {
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                if (visible) {
                    mFilterSwitch.setVisibility(View.VISIBLE);
                    mExportButton.setVisibility(View.VISIBLE);
                    mEegValsTextView.setVisibility(View.VISIBLE);
                } else {
                    mExportButton.setVisibility(View.INVISIBLE);
                    mFilterSwitch.setVisibility(View.INVISIBLE);
                    mEegValsTextView.setVisibility(View.INVISIBLE);
                }
            }
        });
    }

    private int batteryLevel = -1;
    @Override
    public void onCharacteristicRead(BluetoothGatt gatt, BluetoothGattCharacteristic characteristic, int status) {
        Log.i(TAG, "onCharacteristicRead");
        if (status == BluetoothGatt.GATT_SUCCESS) {
            if (AppConstant.CHAR_BATTERY_LEVEL.equals(characteristic.getUuid())) {
                batteryLevel = characteristic.getIntValue(BluetoothGattCharacteristic.FORMAT_UINT8, 0);
                updateBatteryStatus(batteryLevel, batteryLevel + " %");
                Log.i(TAG, "Battery Level :: " + batteryLevel);
            }
        } else {
            Log.e(TAG, "onCharacteristic Read Error" + status);
        }
    }
    private int dataCnt1000 = 0;
    private boolean eeg_ch1_data_on = false;
    private boolean eeg_ch2_data_on = false;
    private boolean eeg_ch3_data_on = false;
    private boolean eeg_ch4_data_on = false;
    //most recent eeg data packet:
    private int[] eeg_ch1_data = new int[6];
    private int[] eeg_ch2_data = new int[6];
    private int[] eeg_ch3_data = new int[6];
    private int[] eeg_ch4_data = new int[6];
    @Override
    public void onCharacteristicChanged(BluetoothGatt gatt, BluetoothGattCharacteristic characteristic) {
        //TODO: ADD BATTERY MEASURE CAPABILITY IN FIRMWARE: (ble_ADC)
        if (AppConstant.CHAR_BATTERY_LEVEL.equals(characteristic.getUuid())) {
            batteryLevel = characteristic.getIntValue(BluetoothGattCharacteristic.FORMAT_UINT8, 0);
            updateBatteryStatus(batteryLevel, batteryLevel + " %");
            String timeStamp = getTimeStamp2();
            exportLogFile(false, "Battery Level Changed at " + timeStamp + getDetails()+"\r\n");
            Log.i(TAG, "Battery Level :: " + batteryLevel);
        }

        if (AppConstant.CHAR_EEG_CH1_SIGNAL.equals(characteristic.getUuid())) {
            byte[] dataEmgBytes = characteristic.getValue();
            if(!eeg_ch1_data_on) {
                eeg_ch1_data_on = true;
            }
            int byteLength = dataEmgBytes.length;
            //TODO: Remember to check/uncheck plotImplicitXVals (boolean)
            getDataRateBytes(byteLength);
            for (int i = 0; i < byteLength/3; i++) { //0→9
                dataCnt1000++; //count?
                int data = unsignedBytesToInt(dataEmgBytes[3*i], dataEmgBytes[3*i+1], dataEmgBytes[3*i+2]);
                eeg_ch1_data[i] = unsignedToSigned(data, 24);
                updateEEG(eeg_ch1_data[i]);
            }
//            Log.e("Ch1 = ",String.valueOf(byteLength)+" # of bytes");
        }

        if (AppConstant.CHAR_EEG_CH2_SIGNAL.equals(characteristic.getUuid())) {
            if(!eeg_ch2_data_on) {
                eeg_ch2_data_on = true;
            }
            byte[] dataEmgBytes = characteristic.getValue();
            int byteLength = dataEmgBytes.length;
            //TODO: Remember to check/uncheck plotImplicitXVals (boolean)
            getDataRateBytes(byteLength);
            for (int i = 0; i < byteLength/3; i++) { //0→9
                dataCnt1000++; //count?
                int data = unsignedBytesToInt(dataEmgBytes[3*i], dataEmgBytes[3*i+1], dataEmgBytes[3*i+2]);
                eeg_ch2_data[i] = unsignedToSigned(data, 24);
            }
//            Log.e("Ch2 = ",String.valueOf(byteLength)+" # of bytes");
        }

        if (AppConstant.CHAR_EEG_CH3_SIGNAL.equals(characteristic.getUuid())) {
            if(!eeg_ch3_data_on) {
                eeg_ch3_data_on = true;
            }
            byte[] dataEmgBytes = characteristic.getValue();
            int byteLength = dataEmgBytes.length;
            //TODO: Remember to check/uncheck plotImplicitXVals (boolean)
            getDataRateBytes(byteLength);
            for (int i = 0; i < byteLength/3; i++) { //0→9
                dataCnt1000++; //count?
                int data = unsignedBytesToInt(dataEmgBytes[3*i], dataEmgBytes[3*i+1], dataEmgBytes[3*i+2]);
                eeg_ch3_data[i] = unsignedToSigned(data, 24);
            }
//            Log.e("Ch3 = ",String.valueOf(byteLength)+" # of bytes");
        }

        if (AppConstant.CHAR_EEG_CH4_SIGNAL.equals(characteristic.getUuid())) {
            if(!eeg_ch4_data_on) {
                eeg_ch4_data_on = true;
            }
            byte[] dataEmgBytes = characteristic.getValue();
            int byteLength = dataEmgBytes.length;
            //TODO: Remember to check/uncheck plotImplicitXVals (boolean)
            getDataRateBytes(byteLength);
            for (int i = 0; i < byteLength/3; i++) { //0→9
                dataCnt1000++; //count?
                int data = unsignedBytesToInt(dataEmgBytes[3*i], dataEmgBytes[3*i+1], dataEmgBytes[3*i+2]);
                eeg_ch4_data[i] = unsignedToSigned(data, 24);
            }
//            Log.e("Ch4 = ",String.valueOf(byteLength)+" # of bytes");
        }
        if(eeg_ch4_data_on && eeg_ch3_data_on && eeg_ch2_data_on && eeg_ch1_data_on) {
            eeg_ch1_data_on = false;
            eeg_ch2_data_on = false;
            eeg_ch3_data_on = false;
            eeg_ch4_data_on = false;
            for (int i = 0; i < 6; i++) {
                writeToDisk24(eeg_ch1_data[i],eeg_ch2_data[i],eeg_ch3_data[i],eeg_ch4_data[i]);
            }
            Log.d(TAG,"Arrays 1-4 Received");
        }
    }


    private int ECGBDSize = 1000;
    private int smallWindowSize = 250;
    private int afibIndex = 0;
    /**
     * TODO: FIX THIS:
     */
    private int unfiltIndex = 0;
    private double[] unfilteredEcgSignal = new double[250]; //250 or 500
    private double[] ECGBufferUnfiltered = new double[ECGBDSize]; //1000
    private double[] ECGBufferFiltered2 = new double[ECGBDSize];
//    private double[] ECGBufferAfibDetection = new double[3000];

    public void sendSMS(String phoneNumber, String msg) {
        try {
            SmsManager smsManager = SmsManager.getDefault();
            smsManager.sendTextMessage(phoneNumber, null, msg, null, null);
            Toast.makeText(getApplicationContext(), "Fall Alert Sent to Healthcare Provider", Toast.LENGTH_LONG).show();
        } catch (Exception ex) {
            Toast.makeText(getApplicationContext(),ex.getMessage(),Toast.LENGTH_LONG).show();
            ex.printStackTrace();
        }
    }

    private void writeToDisk(final int value) {
        //Add to Back:
        double dividedInt = (double) value / 32767.0;
        double dataVoltage = (dividedInt * 2.42);
        //TODO: Exporting unfiltered data to drive.
        try {
            exportFile(false, false, "", dataVoltage);
        } catch (IOException e) {
            Log.e("IOException", e.toString());
        }
    }
    private void writeToDisk24(final int value) {
        //Add to Back:
//        double dividedInt = (double) value / 8388607.0;
//        double dataVoltage = (dividedInt * 2.42);
        //TODO: Exporting unfiltered data to drive.
        try {
            exportFile(false, false, "", convert24bitInt(value));
        } catch (IOException e) {
            Log.e("IOException", e.toString());
        }
    }

    private void writeToDisk24(final int ch1, final int ch2, final int ch3, final int ch4) {
        //Add to Back: TODO: Exporting unfiltered data to drive.
        try {
            for (int i = 0; i < 6; i++) {
                double ch1d = convert24bitInt(ch1);
                double ch2d = convert24bitInt(ch2);
                double ch3d = convert24bitInt(ch3);
                double ch4d = convert24bitInt(ch4);
                exportFile(false, false, "", ch1d,ch2d,ch3d,ch4d);
            }
        } catch (IOException e) {
            Log.e("IOException", e.toString());
        }
    }

    private double convert24bitInt(final int int24bit) {
        double dividedInt = (double) int24bit/8388607.0;
        return dividedInt*2.42;
    }

    private int unfiltIndexCh1 = 0;
    private double[] explicitXValsCh1 = new double[1000];
    private double[] explicitXValsCh2 = new double[1000];
    private double[] explicitXValsCh3 = new double[1000];
    private double[] explicitXValsCh4 = new double[1000];
    private double[] unfiltEEGCH1 = new double[250];
    private double[] EEGBufferCh1Unfiltered = new double[1000];
    private double[] EEGBufferCh1Filtered = new double[1000];
    private int unfiltIndexCh2 = 0;
    private double[] unfiltEEGCH2 = new double[250];
    private double[] EEGBufferCh2Unfiltered = new double[1000];
    private double[] EEGBufferCh2Filtered = new double[1000];
    private int unfiltIndexCh3 = 0;
    private double[] unfiltEEGCH3 = new double[250];
    private double[] EEGBufferCh3Unfiltered = new double[1000];
    private double[] EEGBufferCh3Filtered = new double[1000];
    private int unfiltIndexCh4 = 0;
    private double[] unfiltEEGCH4 = new double[250];
    private double[] EEGBufferCh4Unfiltered = new double[1000];
    private double[] EEGBufferCh4Filtered = new double[1000];

   /* private void updateEEG(final int value, final int channelNum) {
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                double dividedInt = (double) value / 8388607.0;
                double dataVoltage = (dividedInt * 2.42);
                switch (channelNum) {
                    case 1:
                        mEegValsTextView.setText(" " + " IntVal="+String.valueOf(value));
                        unfiltEEGCH1[unfiltIndexCh1] = dataVoltage;
                        if((unfiltIndexCh1 % 249==0) && unfiltIndexCh1!=0) {
                            for (int i = 0; i < 750; i++) {
                                //Shift Back EEGBufferCh1Unfiltered
                                EEGBufferCh1Unfiltered[i] = EEGBufferCh1Unfiltered[i+250];
                                //Shift Back ExplicitX;
                                explicitXValsLong[i] = explicitXValsLong[250+i];
                            }
                            for (int i = 0; i < 250; i++) {
                                //Add to front (ECGBufferUnfiltered):
                                EEGBufferCh1Unfiltered[750+i] = unfiltEEGCH1[i];
                                //Add to front (explicitX)
                                explicitXValsLong[750+i] = ((double) mTimeSeconds) + ((double)i/250);
                            }
                            EEGBufferCh1Filtered = jBwFilter(EEGBufferCh1Unfiltered);
                            if(eegDataSeries1.size()>249) {
                                double max = findGraphMax(eegDataSeries1, eegDataSeries2, eegDataSeries3, eegDataSeries4);
                                double min = findGraphMin(eegDataSeries1, eegDataSeries2, eegDataSeries3, eegDataSeries4);
                                if((max-min)!=0) {
                                    if(currentBM!=BoundaryMode.AUTO) {
                                        eegPlot.setRangeBoundaries(-2.5, 2.5, BoundaryMode.AUTO);
                                        currentBM = BoundaryMode.AUTO;
                                    }
                                    eegPlot.setRangeStepValue((max-min)/5);
                                } else {
                                    if(currentBM!=BoundaryMode.FIXED) {
                                        eegPlot.setRangeBoundaries(min-1, max+1, BoundaryMode.FIXED);
                                        currentBM = BoundaryMode.FIXED;
                                    }
                                    eegPlot.setRangeStepValue(2.0/5.0);
                                }
                            }
                            //TODO: Now plot:

                        } else {
                            unfiltIndexCh1++;
                        }
                        break;
                    case 2:
                        unfiltEEGCH2[unfiltIndexCh2] = dataVoltage;
                        if((unfiltIndexCh2 % 249==0) && unfiltIndexCh2!=0) {
                            //TODO:
                        } else {
                            unfiltIndexCh2++;
                        }
                        break;
                    case 3:
                        unfiltEEGCH3[unfiltIndexCh3] = dataVoltage;
                        if((unfiltIndexCh3 % 249==0) && unfiltIndexCh3!=0) {
                            //TODO:
                        } else {
                            unfiltIndexCh3++;
                        }
                        break;
                    case 4:
                        unfiltEEGCH4[unfiltIndexCh4] = dataVoltage;
                        if((unfiltIndexCh4 % 249==0) && unfiltIndexCh4!=0) {
                            //TODO:
                        } else {
                            unfiltIndexCh4++;
                        }
                        break;
                    default:
                        break;
                }
            }
        });
    }*/

    private void updateEEG(final int value) {
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                //Add to Back:
                mEegValsTextView.setText(" " + " IntVal="+String.valueOf(value));
                double dividedInt = (double) value / 8388607.0;
                double dataVoltage = (dividedInt * 2.42);
                //TODO: Exporting unfiltered data to drive.
                unfilteredEcgSignal[unfiltIndex] = dataVoltage;

                if((unfiltIndex % 249==0) && unfiltIndex!=0) { //SAMPLES_PER_SECOND - 1
                    for (int i = 0; i < 750; i++) {
                        //Shift Back ECGBuffer2
                        ECGBufferUnfiltered[i] = ECGBufferUnfiltered[i+250];
                        //Shift Back ExplicitX;
                        explicitXValsLong[i] = explicitXValsLong[250+i];
                    }

                    for (int i = 0; i < 250; i++) {
                        //Add to front (ECGBufferUnfiltered):
                        ECGBufferUnfiltered[750+i] = unfilteredEcgSignal[i];
                        //Add to front (explicitX)
                        explicitXValsLong[750+i] = ((double) mTimeSeconds) + ((double)i/250);
                    }
                    //filter:
//                    ECGBufferFiltered2 = jBwFilter(ECGBufferUnfiltered);
                    ECGBufferFiltered2 = jeegcfilt(ECGBufferUnfiltered);
                    //Todo: adjust Range Step Value every 1 s:
                    //Every 4 seconds elapsed:
                    //TODO: Now Analyze filtered data (after plot)
                    if(eegDataSeries1.size()>249) {
                        double max = findGraphMax(eegDataSeries1);
                        double min = findGraphMin(eegDataSeries1);
                        if((max-min)!=0) {
                            if(currentBM!=BoundaryMode.AUTO) {
                                eegPlot.setRangeBoundaries(-2.5, 2.5, BoundaryMode.AUTO);
                                currentBM = BoundaryMode.AUTO;
                            }
                            eegPlot.setRangeStepValue((max-min)/5);
                        } else {
                            if(currentBM!=BoundaryMode.FIXED) {
                                eegPlot.setRangeBoundaries(min-1, max+1, BoundaryMode.FIXED);
                                currentBM = BoundaryMode.FIXED;
                            }
                            eegPlot.setRangeStepValue(2.0/5.0);
                        }
                    }
                    //TODO: Now plot:
                    //The idea is that every time this thing is triggered, we want to plot the entirety
                    //of [explicitXValsLong, ECGBufferFiltered2] BEFORE the next second occurs.
                    if(mTimeSeconds == HISTORY_SECONDS) {
                        newMinX = Math.floor(explicitXValsLong[0]);
                        newMaxX = Math.floor(explicitXValsLong[999]);
                        eegPlot.setDomainBoundaries(newMinX, newMaxX, BoundaryMode.AUTO);
                    }
                    /*
                     * TODO: Probably set graph so that it plots at 1ms.??
                     * clearGraph(), and plotGraph() work, but doesn't look good
                     */
                    if(!plotImplicitXVals) {
                        if (!startedGraphingEcgData) {
                            if (!initExecutor) {
                                startGraphEcgData3();
                                initExecutor = true;
                            }
                            startedGraphingEcgData = true;
                        }
                    }
                    mTimeSeconds++;
                    unfiltIndex=0;
                } else {
                    unfiltIndex++;
                }
//                mEegValsTextView.setText(" " + " Voltage="+String.format("%1.4f",dataVoltage)+"V");
            }
        });
    }

    /**
     * ### This method is for <b>250SPS</b> only
     * @param value Single ECG value
     */
    private void shortUpdateECGState(final int value) {
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
            //Add to Back:
            double dividedInt = (double) value / 32767.0;
            double dataVoltage = (dividedInt * 2.42);
            //TODO: Exporting unfiltered data to drive.
            unfilteredEcgSignal[unfiltIndex] = dataVoltage;

            if((unfiltIndex % 249==0) && unfiltIndex!=0) { //SAMPLES_PER_SECOND - 1
                for (int i = 0; i < 750; i++) {
                    //Shift Back ECGBuffer2
                    ECGBufferUnfiltered[i] = ECGBufferUnfiltered[i+250];
                    //Shift Back ExplicitX;
                    explicitXValsLong[i] = explicitXValsLong[250+i];
                }

                for (int i = 0; i < 250; i++) {
                    //Add to front (ECGBufferUnfiltered):
                    ECGBufferUnfiltered[750+i] = unfilteredEcgSignal[i];
                    //Add to front (explicitX)
                    explicitXValsLong[750+i] = ((double) mTimeSeconds) + ((double)i/250);
                }
                //filter:
//                ECGBufferFiltered2 = jFirFilter(ECGBufferUnfiltered);
                ECGBufferFiltered2 = jBwFilter(ECGBufferUnfiltered);
                //Todo: adjust Range Step Value every 1 s:
                //Every 4 seconds elapsed:
                //TODO: Now Analyze filtered data (after plot)
                if(eegDataSeries1.size()>249) {
                    double max = findGraphMax(eegDataSeries1);
                    double min = findGraphMin(eegDataSeries1);
                    if((max-min)!=0) {
                        if(currentBM!=BoundaryMode.AUTO) {
                            eegPlot.setRangeBoundaries(-2.5, 2.5, BoundaryMode.AUTO);
                            currentBM = BoundaryMode.AUTO;
                        }
                        eegPlot.setRangeStepValue((max-min)/5);
                    } else {
                        if(currentBM!=BoundaryMode.FIXED) {
                            eegPlot.setRangeBoundaries(min-1, max+1, BoundaryMode.FIXED);
                            currentBM = BoundaryMode.FIXED;
                        }
                        eegPlot.setRangeStepValue(2.0/5.0);
                    }
                }
                //TODO: Now plot:
                 //The idea is that every time this thing is triggered, we want to plot the entirety
                 //of [explicitXValsLong, ECGBufferFiltered2] BEFORE the next second occurs.
                if(mTimeSeconds == HISTORY_SECONDS) {
                    newMinX = Math.floor(explicitXValsLong[0]);
                    newMaxX = Math.floor(explicitXValsLong[999]);
                    eegPlot.setDomainBoundaries(newMinX, newMaxX, BoundaryMode.AUTO);
                }
                /**
                 * TODO: Probably set graph so that it plots at 1ms.??
                 * clearGraph(), and plotGraph() work, but doesn't look good
                 */
                if(!plotImplicitXVals) {
                    if (!startedGraphingEcgData) {
                        if (!initExecutor) {
                            startGraphEcgData3();
                            initExecutor = true;
                        }
                        startedGraphingEcgData = true;
                    }
                }
                mTimeSeconds++;
                unfiltIndex=0;
            } else {
                unfiltIndex++;
            }
            mEegValsTextView.setText(" " + " Voltage="+String.format("%1.4f",dataVoltage)+"V");
            }
        });
    }

    private double findGraphMax(SimpleXYSeries s1, SimpleXYSeries s2, SimpleXYSeries s3, SimpleXYSeries s4) {
        double max1 = (double) s1.getY(0);
        double max2 = (double) s2.getY(0);
        double max3 = (double) s3.getY(0);
        double max4 = (double) s4.getY(0);
        //TODO: smallest array size???
        for (int i = 1; i < s1.size(); i++) {
            double a1 = (double) s1.getY(i);
            if (a1 > max1) {
                max1 = a1;
            }
        }
        for (int i = 0; i < s2.size(); i++) {
            double a2 = (double) s2.getY(i);
            if (a2>max2) {
                max2 = a2;
            }
        }
        for (int i = 0; i < s3.size(); i++) {
            double a3 = (double) s3.getY(i);
            if (a3>max3) {
                max3 = a3;
            }
        }
        for (int i = 0; i < s4.size(); i++) {
            double a4 = (double) s4.getY(i);
            if (a4>max4) {
                max4 = a4;
            }
        }
        double[] maxVals = {max1, max2, max3, max4};
        Arrays.sort(maxVals);
        return maxVals[1];
    }

    private double findGraphMax(SimpleXYSeries s) {
        double max = (double)s.getY(0);
        for (int i = 1; i < s.size(); i++) {
            double a = (double)s.getY(i);
            if(a>max) {
                max = a;
            }
        }
        return max;
    }

    private double findGraphMin (SimpleXYSeries s1, SimpleXYSeries s2, SimpleXYSeries s3, SimpleXYSeries s4) {
        double min1 = (double) s1.getY(0);
        double min2 = (double) s2.getY(0);
        double min3 = (double) s3.getY(0);
        double min4 = (double) s4.getY(0);
        for (int i = 1; i < s1.size(); i++) {
            double a1 = (double) s1.getY(i);
            if (a1 < min1) {
                min1 = a1;
            }
        }
        for (int i = 0; i < s2.size(); i++) {
            double a2 = (double) s2.getY(i);
            if (a2 < min2) {
                min2 = a2;
            }
        }
        for (int i = 0; i < s3.size(); i++) {
            double a3 = (double) s3.getY(i);
            if (a3 < min3) {
                min3 = a3;
            }
        }
        for (int i = 0; i < s4.size(); i++) {
            double a4 = (double) s4.getY(i);
            if (a4 < min4) {
                min4 = a4;
            }
        }
        double[] minVals = {min1,min2, min3, min4};
        Arrays.sort(minVals);
        return minVals[4];
    }

    private double findGraphMin(SimpleXYSeries s) {
        double min = (double)s.getY(0);
        for (int i = 1; i < s.size(); i++) {
            double a = (double)s.getY(i);
            if(a<min) {
                min = a;
            }
        }
        return min;
    }



    private void stopGraphEcgData() {
        startedGraphingEcgData = false;
    }

    private void startGraphEcgData3() {
        executor.scheduleAtFixedRate(new graphEcgData3(), 0, periodShort, TimeUnit.MICROSECONDS);
        startedGraphingEcgData = true;
    }

    private class graphEcgData3 implements Runnable {
        @Override
        public void run() {
            if (startedGraphingEcgData) {
                if (eegIndex < 1000) {
                    if (eegDataSeries1.size() > 1000) {
                        eegDataSeries1.removeFirst();
                    }
                    if (plotImplicitXVals) {
                        if (filterData) {
                            eegDataSeries1.addLast(null, ECGBufferFiltered2[eegIndex]);
                            /*if(!filterType){
                            } else {
                                eegDataSeries1.addLast(null, ECGBufferFiltered3[eegIndex]);
                            }*/
                        } else {
                            eegDataSeries1.addLast(null, ECGBufferUnfiltered[eegIndex]);
                        }
                    } else {
                        if (filterData) {
                            eegDataSeries1.addLast(explicitXValsLong[eegIndex], ECGBufferFiltered2[eegIndex]);
                            /*if(!filterType){
                            } else {
                                eegDataSeries1.addLast(explicitXValsLong[eegIndex], ECGBufferFiltered3[eegIndex]);
                            }*/
                        } else {
                            eegDataSeries1.addLast(explicitXValsLong[eegIndex], ECGBufferUnfiltered[eegIndex]);
                        }
                    }
                    eegIndex++;
                } else {
                    eegIndex = 750;
                    lastFilteredPacket++;
                    stopGraphEcgData();
                    Log.d("Graphing Packet#", String.valueOf(lastFilteredPacket));
                }
            }
        }
    }

    private double[] ECGBufferDoubles = new double[ECGBDSize]; //4s of data, unfiltered
    private double[] ECGBufferFiltered = new double[ECGBDSize];
    private double[] explicitXValsLong = new double[ECGBDSize];
    private boolean initExecutor = false;
    private boolean startedGraphingEcgData = false;

    private long lastFilteredPacket = 0;

    //TODO: Display results at fixed rate:

    private Number newMinX;
    private Number newMaxX;

    private String getDetails() {
        return "Details:\r\n" +
                "Bluetooth LE Device Name: "+mBluetoothDevice.getName() + "\r\n" +
                "Bluetooth LE Address: "+mBluetoothDevice.getAddress() + "\r\n" +
                "Last Battery Level: "+String.valueOf(batteryLevel)+"% \r\n\r\n" +
                "Last RSSI: "+mLastRssi + "\r\n" +
                "Sampling Rate: "+String.valueOf(DATA_RATE_SAMPLES_PER_SECOND)+"Hz\r\n" +
                "Samples Read: "+String.valueOf(dataCnt1000)+ "\r\n" +
                "Android Device Manufacturer: "+ Build.MANUFACTURER + "\r\n" +
                "Device Model:" + Build.MODEL + "\r\n" +
                "Android Unique ID: " + Settings.Secure.getString(getBaseContext()
                .getContentResolver(), Settings.Secure.ANDROID_ID) + "\r\n" +
                "SDK Version: " + String.valueOf(Build.VERSION.SDK_INT)+ "\r\n" +
                "Version Release: " + Build.VERSION.RELEASE + "\r\n" +
                "Android Battery Level: " + androidDeviceBatteryLevel + "\r\n" +
                "Battery Status: " + androidDeviceBatteryStatus + "\r\n" + "\r\n" +
                "";
    }
    private void getDataRateBytes(int bytes) {
        mCurrentTime = System.currentTimeMillis();
        points += bytes;
        if (mCurrentTime > (mLastTime + 5000)) {
            dataRate = (points / 5);
            points = 0;
            mLastTime = mCurrentTime;
            Log.e(" DataRate:", String.valueOf(dataRate) + " Bytes/s");
            runOnUiThread(new Runnable() {
                @Override
                public void run() {

                    mDataRate.setText(String.valueOf(dataRate)+ " Bytes/s");
                }
            });
        }
    }
    //Get Data Rate (assigned to ECG)
    private void getDataRate(int pointnum) {
        mCurrentTime = System.currentTimeMillis();
        points += pointnum;
        if (mCurrentTime > (mLastTime + 5000)) {
            dataRate = (points / 5);
            points = 0;
            mLastTime = mCurrentTime;
            Log.e("ECG DataRate:", String.valueOf(dataRate) + " Hz");
            runOnUiThread(new Runnable() {
                @Override
                public void run() {

                    mDataRate.setText(String.valueOf(dataRate)+ " Hz");
                }
            });
        }
    }

    //Get Data Rate 2 (assigned to MPU)
    private void getDataRate2() {
        mCurrentTime2 = System.currentTimeMillis();
        dataCount++;
        if (mCurrentTime2 > (mLastTime2 + 5000)) {
            dataRate2 = (dataCount / 5);
            dataCount = 0;
            mLastTime2 = mCurrentTime2;
            Log.e("MPU Data Rate:", String.valueOf(dataRate2) + " Hz");
        }
    }

    private int[] getIntArrayFromByteArrayEcg(byte[] rawValue) {
        int[] ecgArray = new int[rawValue.length/2];
        for (int i = 0; i < (rawValue.length/2); i++) {
            byte[] bin = {rawValue[2 * i], rawValue[2 * i + 1]};
//            byte[] bin = {rawValue[2 * i + 1], rawValue[2 * i]};
            ecgArray[i] = getIntfromByte(bin);
        }
        return ecgArray;
    }

    /**
     * @param rawValue - Received from BLE Characteristic
     * @return mpuArray - signed integer values
     */
    private int[] getIntArrFromByteArr(byte[] rawValue) {
        int size = rawValue.length/2;
        int[] mpuArray = new int[6];
        for (int i = 0; i < size; i++) {
            byte[] bin = {rawValue[2 * i], rawValue[2 * i + 1]};
            mpuArray[i] = getIntfromByte(bin);
        }
        return mpuArray;
    }

    /**
     * This function receives a signed 16-bit integer (2-bytes).
     *
     * @param rawValue (byte array)
     * @return -> Signed *int* value
     */
    private int getIntfromByte(byte[] rawValue) {//Assumes 2 bytes
        return (int) rawValue[0] + ((int) rawValue[1] << 8);
    }

    /**
     * Ignore this, it basically does the same as the above funcs.
     * @param rawValue byte array of type int
     * @return intValue returns a single integer value up to 4-bytes (32-bit int)
     */
    private int getIntfromByteAlt(byte[] rawValue) {//Up to 4 bytes
        int intValue = 0;
        if (rawValue.length > 0) intValue = (int) rawValue[0];
        if (rawValue.length > 1) intValue = intValue + ((int) rawValue[1] << 8); //16-bit
        if (rawValue.length > 2) intValue = intValue + ((int) rawValue[2] << 8);
        if (rawValue.length > 3) intValue = intValue + ((int) rawValue[3] << 8);
        return intValue;
    }

    /**
     * Convert an unsigned integer value to a two's-complement encoded
     * signed value.
     */
    private int unsignedToSigned(int unsigned, int size) {
        if ((unsigned & (1 << size - 1)) != 0) {
            unsigned = -1 * ((1 << size - 1) - (unsigned & ((1 << size - 1) - 1)));
        }
        return unsigned;
    }

    /**
     * Convert a signed byte to an unsigned int.
     */
    private int unsignedByteToInt(byte b) {
        return b & 0xFF;
    }

    /**
     * Convert signed bytes to a 16-bit unsigned int.
     */
    private int unsignedBytesToInt(byte b0, byte b1) {
        return (unsignedByteToInt(b0) + (unsignedByteToInt(b1) << 8));
    }

    private int unsignedBytesToInt(byte b0, byte b1, byte b2) {
        return (unsignedByteToInt(b0) + (unsignedByteToInt(b1) << 8) + (unsignedByteToInt(b2) << 16));
    }


    private static final char[] HEX_CHARS = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A',
            'B', 'C', 'D', 'E', 'F' };
    public static String toHexString(byte[] bytes) {
        char[] hexChars = new char[bytes.length * 2];
        int v;
        for (int j = 0; j < bytes.length; j++) {
            v = bytes[j] & 0xFF;
            hexChars[j * 2] = HEX_CHARS[v >>> 4];
            hexChars[j * 2 + 1] = HEX_CHARS[v & 0x0F];
        }
        return new String(hexChars);
    }

    @Override
    public void onReadRemoteRssi(BluetoothGatt gatt, int rssi, int status) {
        uiRssiUpdate(rssi);
        mLastRssi = String.valueOf(rssi)+"db";
    }

    @Override
    public void onConnectionStateChange(BluetoothGatt gatt, int status, int newState) {
        switch (newState) {
            case BluetoothProfile.STATE_CONNECTED:
                mConnected = true;
                runOnUiThread(new Runnable() {
                    @Override
                    public void run() {
                        if(menu!=null) {
                            menu.findItem(R.id.menu_connect).setVisible(false);
                            menu.findItem(R.id.menu_disconnect).setVisible(true);
                        }
                    }
                });
                Log.i(TAG, "Connected");
                updateConnectionState(getString(R.string.connected));
                invalidateOptionsMenu();
                runOnUiThread(new Runnable() {
                    @Override
                    public void run() {
                        mDataRate.setTextColor(Color.BLACK);
                        mDataRate.setTypeface(null, Typeface.NORMAL);
                    }
                });
                exportLogFile(false, "Connected @ "+getTimeStamp()+"\r\n"+"\r\n");
                //Start the service discovery:
                gatt.discoverServices();
                startMonitoringRssiValue();
                break;
            case BluetoothProfile.STATE_DISCONNECTED:
                mConnected = false;
                runOnUiThread(new Runnable() {
                    @Override
                    public void run() {
                        if(menu!=null) {
                            menu.findItem(R.id.menu_connect).setVisible(true);
                            menu.findItem(R.id.menu_disconnect).setVisible(false);
                        }
                    }
                });
                Log.i(TAG, "Disconnected");
                runOnUiThread(new Runnable() {
                    @Override
                    public void run() {
                        mDataRate.setTextColor(Color.RED);
                        mDataRate.setTypeface(null, Typeface.BOLD);
                        mDataRate.setText("0 Hz");
                    }
                });
                exportLogFile(false, "Disconnected @ "+getTimeStamp()+"\r\n"+"\r\n");
                //TODO: ATTEMPT TO RECONNECT:
                updateConnectionState(getString(R.string.disconnected));
                stopMonitoringRssiValue();
                stopGraphEcgData();
                invalidateOptionsMenu();
                break;
            default:
                break;
        }
    }

    public void startMonitoringRssiValue() {
        readPeriodicallyRssiValue(true);
    }

    public void stopMonitoringRssiValue() {
        readPeriodicallyRssiValue(false);
    }

    public void readPeriodicallyRssiValue(final boolean repeat) {
        mTimerEnabled = repeat;
        // check if we should stop checking RSSI value
        if (!mConnected || mBluetoothGatt == null || !mTimerEnabled) {
            mTimerEnabled = false;
            return;
        }

        mTimerHandler.postDelayed(new Runnable() {
            @Override
            public void run() {
                if (mBluetoothGatt == null  || !mConnected) {
                    mTimerEnabled = false;
                    return;
                }

                // request RSSI value
                mBluetoothGatt.readRemoteRssi();
                // add call it once more in the future
                readPeriodicallyRssiValue(mTimerEnabled);
            }
        }, RSSI_UPDATE_TIME_INTERVAL);
    }

    @Override
    public void onCharacteristicWrite(BluetoothGatt gatt, BluetoothGattCharacteristic
            characteristic, int status) {
        Log.i(TAG, "onCharacteristicWrite :: Status:: " + status);
    }

    @Override
    public void onDescriptorWrite(BluetoothGatt gatt, BluetoothGattDescriptor descriptor, int status) {
    }

    @Override
    public void onDescriptorRead(BluetoothGatt gatt, BluetoothGattDescriptor descriptor, int status) {
        Log.i(TAG, "onDescriptorRead :: Status:: " + status);
    }

    @Override
    public void onError(String errorMessage) {
        Log.e(TAG, "Error:: " + errorMessage);
    }

    private void updateConnectionState(final String status) {
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                if (status.equals(getString(R.string.connected))) {
                    Toast.makeText(getApplicationContext(), "Device Connected!", Toast.LENGTH_SHORT).show();
                } else if (status.equals(getString(R.string.disconnected))) {
                    Toast.makeText(getApplicationContext(), "Device Disconnected!", Toast.LENGTH_SHORT).show();
                }
            }
        });
    }

    private void batteryNotAvailable() {
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
//                Toast.makeText(getApplicationContext(), "Device Does Not Have Battery", Toast.LENGTH_SHORT).show();
            }
        });
    }

    private double[] explicitXVals = new double[250];
    private double[] responseFilteredEcg = new double[250];

    private void updateBatteryStatus(final int percent, final String status) {
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                if (percent <= batteryWarning) {
                    mBatteryLevel.setTextColor(Color.RED);
                    mBatteryLevel.setTypeface(null, Typeface.BOLD);
                    Toast.makeText(getApplicationContext(), "Charge Battery, Battery Low " + status, Toast.LENGTH_SHORT).show();
                } else {
                    mBatteryLevel.setTextColor(Color.GREEN);
                    mBatteryLevel.setTypeface(null, Typeface.BOLD);
                }
                mBatteryLevel.setText(status);
            }
        });
    }

    private void uiRssiUpdate(final int rssi) {
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                MenuItem menuItem = menu.findItem(R.id.action_rssi);
                MenuItem status_action_item = menu.findItem(R.id.action_status);
                final String valueOfRSSI = String.valueOf(rssi) + " dB";
//                mRssi.setText(valueOfRSSI);
                menuItem.setTitle(valueOfRSSI);
                if (mConnected) {
                    String newStatus = "Status: " + getString(R.string.connected);
                    status_action_item.setTitle(newStatus);
                } else {
                    String newStatus = "Status: " + getString(R.string.disconnected);
                    status_action_item.setTitle(newStatus);
                }
            }
        });
    }

    /*
    * Application of JNI code:
    */
    static {
        System.loadLibrary("android-jni");
    }

    //ECG FIR FILTER:
    private native int jmainFHC(boolean b);

    private native double[] jfullHybridClassifier(double[] data1, double[] data2, double[] data3, double[] data4, boolean EOGOnly); //size = 1000

    //ECG BW Filter:
    private native int jmainEegFilt(boolean b);

    private native double[] jeegcfilt(double[] array);

    private native double[] jBwFilter(double[] ecg);
}
