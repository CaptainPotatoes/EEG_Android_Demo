package com.mahmoodms.bluetooth.ecgfallsensordemo;

import android.app.ActionBar;
import android.app.Activity;
import android.app.NotificationManager;
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
import android.support.v4.app.NotificationCompat;
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
import com.androidplot.xy.BarFormatter;
import com.androidplot.xy.BarRenderer;
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
    private boolean mConnected;
    //Class instance variable
    private BluetoothLe mBluetoothLe;
    private BluetoothManager mBluetoothManager = null;
    private BluetoothGatt mBluetoothGatt = null;
    private BluetoothDevice mBluetoothDevice;
    //Layout - TextViews and Buttons

    private TextView mEcgValues;
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

    private boolean filterData = true;
    private int dataCount = 0;
    private int points = 0;
    private Menu menu;
    //RSSI Stuff:
    private static final int RSSI_UPDATE_TIME_INTERVAL = 2000;
    private static final String BRADYCARDIA_WARNING = "Bradycardia Detected! - Low Heart Rate";
    private static final String TACHYCARDIA_WARNING = "Tachycardia Detected! - High Heart Rate";
    private static final String AFIB_WARNING = "Atrial Fibrillation Detected!";
    private static final String AFLUTTER_WARNING = "Atrial Flutter Detected!";
    private Handler mTimerHandler = new Handler();

    private boolean mTimerEnabled = false;
    /**
     * Initialize Plot:
     */
    private XYPlot ecgPlot;
    private XYPlot accelerometerPlot;
    private Redrawer redrawer;
    private SimpleXYSeries accelerometerDataSeries;
    private SimpleXYSeries plotAccDataSeries;
    private SimpleXYSeries ecgDataSeries;
    private static final int HISTORY_SIZE = 1500;
    private static final int HISTORY_SIZE_2 = 1000;
    private static final int HISTORY_SECONDS = 6;
    private static final int HISTORY_SECONDS_2 = 4;
    private boolean plotImplicitXVals = false;
    private int DATA_RATE_SAMPLES_PER_SECOND = 0;
    //for bounds:
    private BoundaryMode currentBM = BoundaryMode.AUTO;
    //Data Variables:
    private int batteryWarning = 20;//%
    private double[] allMpuData = new double[600]; // First 600 = thetanull.
    private double[] yfit = new double[5];
    private double accelerometer_resultant;
    private int[] mpuData = new int[6];
    final private boolean terminate = true;
    final private boolean initialize = false;
    private String fileTimeStamp = "";
    private long periodShort = 3500; // 3.5ms (changed from 3.7ms on 10/11 @  7:50pm)
    private int ecgIndex = 0;
    private double dataRate;
    private double dataRate2;
    /* Notification stuff */

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
        //Set up action bar:
        if (getActionBar() != null) {
            getActionBar().setDisplayHomeAsUpEnabled(true);
        }
        ActionBar actionBar = getActionBar();
        actionBar.setBackgroundDrawable(new ColorDrawable(Color.parseColor("#6078ef")));

        //Flag to keep screen on (stay-awake):
        getWindow().addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);
        //Set up TextViews
        mEcgValues = (TextView) findViewById(R.id.ecgValue);
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
        accelerometerDataSeries = new SimpleXYSeries("Motion Sensor Data (g)");
        plotAccDataSeries = new SimpleXYSeries("MotionSensorData");
        ecgDataSeries = new SimpleXYSeries("ECG Data (V)");
        ecgPlot = (XYPlot) findViewById(R.id.ecgPlot);
        accelerometerPlot = (XYPlot) findViewById(R.id.accelLevelsPlot);
        //Todo: Graph temporarily uses data index - find alternative to implicit XVals→(seconds)
        if (plotImplicitXVals) {
            ecgDataSeries.useImplicitXVals();
            ecgPlot.setDomainBoundaries(0, HISTORY_SIZE_2, BoundaryMode.FIXED);
            ecgPlot.setDomainStepMode(XYStepMode.INCREMENT_BY_VAL);
            ecgPlot.setDomainStepValue(HISTORY_SIZE_2/5);
        } else {
            ecgPlot.setDomainBoundaries(0, HISTORY_SECONDS_2, BoundaryMode.FIXED);
            ecgPlot.setDomainStepMode(XYStepMode.INCREMENT_BY_VAL);
            ecgPlot.setDomainStepValue(HISTORY_SECONDS_2 / 4);
        }

        //TODO: Adaptive Range?, or some method of figuring out what is connected
        if(filterData) {
            ecgPlot.setRangeBoundaries(-2.5, 2.5, BoundaryMode.AUTO); //EMG only!
            ecgPlot.setRangeStepValue(1);
        } /*else {
            ecgPlot.setRangeBoundaries(-2.5, 2.5, BoundaryMode.FIXED);
            ecgPlot.setRangeStepValue(0.5);
        }*/
        ecgPlot.setRangeStepMode(XYStepMode.INCREMENT_BY_VAL);
        ecgPlot.setDomainLabel("Time (seconds)");
        ecgPlot.getDomainLabelWidget().pack();
        ecgPlot.setRangeLabel("Voltage (mV)");
        ecgPlot.getRangeLabelWidget().pack();
        ecgPlot.setRangeValueFormat(new DecimalFormat("#.###"));
        ecgPlot.setDomainValueFormat(new DecimalFormat("#"));
        ecgPlot.getDomainLabelWidget().getLabelPaint().setColor(Color.BLACK);
        ecgPlot.getDomainLabelWidget().getLabelPaint().setTextSize(20);
        ecgPlot.getRangeLabelWidget().getLabelPaint().setColor(Color.BLACK);
        ecgPlot.getRangeLabelWidget().getLabelPaint().setTextSize(20);
        ecgPlot.getGraphWidget().getDomainTickLabelPaint().setColor(Color.BLACK);
        ecgPlot.getGraphWidget().getRangeTickLabelPaint().setColor(Color.BLACK);
        ecgPlot.getGraphWidget().getDomainTickLabelPaint().setTextSize(23); //TODO: was 36
        ecgPlot.getGraphWidget().getRangeTickLabelPaint().setTextSize(23);
        ecgPlot.getGraphWidget().getDomainGridLinePaint().setColor(Color.WHITE);
        ecgPlot.getGraphWidget().getRangeGridLinePaint().setColor(Color.WHITE);
        ecgPlot.getLegendWidget().getTextPaint().setColor(Color.BLACK);
        ecgPlot.getLegendWidget().getTextPaint().setTextSize(20);
        ecgPlot.getTitleWidget().getLabelPaint().setTextSize(20);
        ecgPlot.getTitleWidget().getLabelPaint().setColor(Color.BLACK);
        LineAndPointFormatter lineAndPointFormatter2 = new LineAndPointFormatter(Color.RED, null, null, null);
        lineAndPointFormatter2.getLinePaint().setStrokeWidth(4);
        LineAndPointFormatter lineAndPointFormatter1 = new LineAndPointFormatter(Color.BLACK, null, null, null);
        lineAndPointFormatter1.getLinePaint().setStrokeWidth(3);
        ecgPlot.addSeries(ecgDataSeries, lineAndPointFormatter1);
        if (plotImplicitXVals) ecgPlot.addSeries(plotAccDataSeries, lineAndPointFormatter2);
        accelerometerPlot.setDomainStepValue(1);
        accelerometerPlot.setTicksPerRangeLabel(1);
        accelerometerPlot.setDomainBoundaries(-1, 1, BoundaryMode.FIXED);
        accelerometerPlot.setRangeBoundaries(0, 3, BoundaryMode.FIXED);
        accelerometerPlot.getGraphWidget().getDomainTickLabelPaint().setColor(Color.TRANSPARENT);
        accelerometerPlot.getGraphWidget().getRangeTickLabelPaint().setColor(Color.BLACK);
        accelerometerPlot.getGraphWidget().getRangeTickLabelPaint().setTextSize(36);
        accelerometerPlot.setDomainLabel("");
        accelerometerPlot.getDomainLabelWidget().pack();
        accelerometerPlot.setRangeLabel("Resultant (g)");
        accelerometerPlot.getRangeLabelWidget().pack();
        accelerometerPlot.setRangeValueFormat(new DecimalFormat("#.#"));
        accelerometerPlot.getRangeLabelWidget().getLabelPaint().setColor(Color.BLACK);
        accelerometerPlot.getDomainLabelWidget().getLabelPaint().setColor(Color.BLACK);
        accelerometerPlot.getRangeLabelWidget().getLabelPaint().setTextSize(20);
        accelerometerPlot.getLegendWidget().getTextPaint().setColor(Color.BLACK);
        accelerometerPlot.getLegendWidget().getTextPaint().setTextSize(20);
        accelerometerPlot.getTitleWidget().getLabelPaint().setTextSize(20);
        accelerometerPlot.getTitleWidget().getLabelPaint().setColor(Color.BLACK);
        accelerometerPlot.addSeries(accelerometerDataSeries,
                new BarFormatter(Color.rgb(200, 0, 0), Color.rgb(0, 80, 0)));

        BarRenderer barRenderer = (BarRenderer) accelerometerPlot.getRenderer(BarRenderer.class);
        if (barRenderer != null) {
            barRenderer.setBarWidth(25);
        }
        redrawer = new Redrawer(
                Arrays.asList(new Plot[]{ecgPlot, accelerometerPlot}),
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
        ecgIndex = 750;
        makeFilterSwitchVisible(false);
        makeFallSensorVisible(false);
        mFilterSwitch.setOnCheckedChangeListener(new CompoundButton.OnCheckedChangeListener() {
            @Override
            public void onCheckedChanged(CompoundButton buttonView, boolean isChecked) {
                /*if(isChecked) {
                    filterData = true;
                } else {
                    filterData = false;
                }*/
                //Simplify using a ternary operator:
//                filterData = (isChecked) ? true : false;
                //Nevermind, I'm an idiot:
                filterData = isChecked;
                if(filterData) {
                    if(ecgDataSeries.size()>0) {
                        double average = 0;
                        for (int i = 0; i < ecgDataSeries.size(); i++) {
                            average+=(double)ecgDataSeries.getY(i);
                            average/=(ecgDataSeries.size());
                        }
                        ecgPlot.setRangeBoundaries(average-0.010, average+0.010, BoundaryMode.AUTO);
                        ecgPlot.setRangeStepValue(0.4);
                    }
                } else {
                    //TODO: ECG is in range of 1-5 mV
                    if(ecgDataSeries.size()>0) {
                        double average = 0;
                        if(ecgDataSeries.size()>250) {
                            for (int i = 0; i < 250; i++) {
                                average+=(double)ecgDataSeries.getY(i);
                                average/=(250.0);
                            }
                        } else {
                            for (int i = 0; i < ecgDataSeries.size(); i++) {
                                average+=(double)ecgDataSeries.getY(i);
                                average/=(ecgDataSeries.size());
                            }
                        }
                        ecgPlot.setRangeBoundaries(average-0.010, average+0.010, BoundaryMode.AUTO);
                        ecgPlot.setRangeStepValue(0.1);
                    }
                }
            }
        });
        mLastTime = System.currentTimeMillis();
        mLastTime2 = mLastTime;
        Arrays.fill(allMpuData, 0.0);
        Arrays.fill(yfit, 0.0);
        Arrays.fill(ECGBufferUnfiltered, 0.0);
        int responsejmainFir = jmainFirFilter(initialize);
        Log.d("ResponseJmainInit:", String.valueOf(responsejmainFir));
        //TODO: Notifications:
        /*Context context = getApplicationContext();
        NotificationCompat.Builder mBuilder = new NotificationCompat.Builder(context)
                .setSmallIcon(R.drawable.ic_launcher)
                .setContentTitle("ECG Notification")
                .setContentText("TestNotification");
        NotificationManager mNotificationManager = (NotificationManager) getSystemService(Context.NOTIFICATION_SERVICE);
        mNotificationManager.notify(mId, mBuilder.build());*/
        this.registerReceiver(this.mBatInfoReceiver, new IntentFilter(Intent.ACTION_BATTERY_CHANGED));
        if(mDeviceName.equals("ECG 1000Hz")) {
            DATA_RATE_SAMPLES_PER_SECOND = 1000;
        } else if (mDeviceName.equals("ECG 500Hz")) {
            DATA_RATE_SAMPLES_PER_SECOND = 500;
        } else if (mDeviceName.equals("ECG 250Hz")||mDeviceName.equals("ECGSensor")) {
            DATA_RATE_SAMPLES_PER_SECOND = 250;
        } else {
            DATA_RATE_SAMPLES_PER_SECOND = 250;
        }
//        mDataRate.setText("Data Rate: "+String.valueOf(DATA_RATE_SAMPLES_PER_SECOND)+"Hz");
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

    //TODO: Write Log File:
    private boolean fileLogInitialized = false;
    private File logFile;
    public void exportLogFile(boolean init, String dataToWrite) {
        if(init) {
            Log.i("exportLogFile", "generated Log file");
            root = Environment.getExternalStorageDirectory();
            File dir = new File(root+"/ECGDataLogs");
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

    private int exportFileDataPointCounter = 0;
    private int exportFilePart = 1;
    public void exportFile(boolean init, boolean terminateExport,
                           String fileName, double ecgData) throws IOException {
        if (init) {
            root = Environment.getExternalStorageDirectory();
            fileTimeStamp = fileName;
            fileExportInitialized = true;
        } else {
            if (fileTimeStamp == null || fileTimeStamp.equals("") || !fileExportInitialized) {
                fileTimeStamp = "ECGSensorData_" + getTimeStamp();
            }
        }
        if (root.canWrite() && init) {
            File dir = new File(root.getAbsolutePath() + "/DataDirectory");
            boolean mkdirsA = dir.mkdirs();
            file = new File(dir, fileTimeStamp + "_part"+ String.valueOf(exportFilePart) + ".csv");
            csvWriter = new CSVWriter(new FileWriter(file));
            Log.d("New File Generated", fileTimeStamp + "_part"+ String.valueOf(exportFilePart) + ".csv");
            exportLogFile(false, "NEW FILE GENERATED: "+fileTimeStamp + "_part"+ String.valueOf(exportFilePart) + ".csv\r\n\r\n");
            if(exportFilePart!=1)exportLogFile(false, getDetails());
        }
        //Write Data to File (if init & terminateExport are both false)
        if (!init && !terminateExport) {
            //TODO: CHANGE THIS TO 1048576
            if(exportFileDataPointCounter<1048575/*15000*/) {
                valueCsvWrite[0] = ecgData + "";
                csvWriter.writeNext(valueCsvWrite);
                exportFileDataPointCounter++;
            } else {
                valueCsvWrite[0] = ecgData + "";
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

    @Override
    public void onResume() {
        makeFilterSwitchVisible(true);
        makeFallSensorVisible(false);
        int responsejmainFir = jmainFirFilter(initialize);
        Log.d("ResponseJmainInit:", String.valueOf(responsejmainFir));
        String fileTimeStampConcat = "ECGSensorData_" + getTimeStamp();
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
    //TODO - remove disconnect from onPause()
    @Override
    protected void onPause() {
        redrawer.pause();
        makeFilterSwitchVisible(false);
        int responsejmainFir = jmainFirFilter(terminate);
        Log.d("ResponseJmainInit:", String.valueOf(responsejmainFir));
        super.onPause();
//        int tempvar = jmainAfibDetectionInit(terminate);
//        stopMonitoringRssiValue();
//        mBluetoothLe.disconnect(mBluetoothGatt);
    }

    private void initializeBluetooth() {
        mBluetoothManager = (BluetoothManager) getSystemService(Context.BLUETOOTH_SERVICE);
        mBluetoothDevice = mBluetoothManager.getAdapter().getRemoteDevice(mDeviceAddress);
        mBluetoothLe = new BluetoothLe(this, mBluetoothManager, this);
        mBluetoothGatt = mBluetoothLe.connect(mBluetoothDevice, false);
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
                connect();
                return true;
            case R.id.menu_disconnect:
                if (mBluetoothLe != null) {
                    if (mBluetoothGatt != null) {
                        mBluetoothGatt.disconnect();
                    }
                }
                return true;
            case android.R.id.home:
                if (mBluetoothLe != null)
                    mBluetoothLe.disconnect(mBluetoothGatt);
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
//                mConnectionState.setText("Connecting...");
                menuItem.setTitle("Connecting...");
//                invalidateOptionsMenu();
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


                if (AppConstant.SERVICE_ION_NA_SIGNAL.equals(service.getUuid())) {
                    mBluetoothLe.setCharacteristicNotification(gatt, service.getCharacteristic(AppConstant.CHAR_ION_NA_SIGNAL), true);
                    batteryNotAvailable();
                }

                if (AppConstant.SERVICE_EMG_SIGNAL.equals(service.getUuid())) {
                    //Set notification for EMG signal:
                    makeFilterSwitchVisible(true);
//                    mBluetoothLe.readCharacteristic(gatt, service.getCharacteristic(AppConstant.CHAR_EMG_DATARATE));
                    mBluetoothLe.setCharacteristicNotification(gatt, service.getCharacteristic(AppConstant.CHAR_EMG_SIGNAL), true);
                    /*for (BluetoothGattCharacteristic characteristic:service.getCharacteristics()) {
                        //TODO: If this characteristic 0x3262
                        if(AppConstant.CHAR_EMG_DATARATE.equals(characteristic.getUuid())) {
                            mBluetoothLe.readCharacteristic(gatt, service.getCharacteristic(AppConstant.CHAR_EMG_DATARATE));
                        }
                    }*/
                }
                if (AppConstant.SERVICE_BATTERY_LEVEL.equals(service.getUuid())) {
                    //Read the device battery percentage
                    mBluetoothLe.readCharacteristic(gatt, service.getCharacteristic(AppConstant.CHAR_BATTERY_LEVEL));
                    mBluetoothLe.setCharacteristicNotification(gatt, service.getCharacteristic(AppConstant.CHAR_BATTERY_LEVEL), true);
                }
                if (AppConstant.SERVICE_MPU.equals(service.getUuid())) {
                    makeFallSensorVisible(true);
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
//                    mEcgAnalysis.setVisibility(View.VISIBLE);
                    mEcgValues.setVisibility(View.VISIBLE);
                } else {
                    mExportButton.setVisibility(View.INVISIBLE);
                    mFilterSwitch.setVisibility(View.INVISIBLE);
//                    mEcgAnalysis.setVisibility(View.INVISIBLE);
                    mEcgValues.setVisibility(View.INVISIBLE);
                }
            }
        });
    }

    private void makeFallSensorVisible(final boolean visible) {
        runOnUiThread(new Runnable() {
            @Override
            public void run() {

                if (visible) {
                    accelerometerPlot.setVisibility(View.VISIBLE);
                } else {
                    accelerometerPlot.setVisibility(View.GONE);
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
            if(AppConstant.CHAR_EMG_DATARATE.equals(characteristic.getUuid())) {
                final byte dataRateCode[] = characteristic.getValue();
                /*if(dataRateCode[0]==(byte)0x01) {
                    DATA_RATE_SAMPLES_PER_SECOND = 250;
                    Log.e("Data Rate:", "250Hz");
                } else if (dataRateCode[0]==(byte)0x02) {
                    DATA_RATE_SAMPLES_PER_SECOND = 500;
                    Log.e("Data Rate:","500Hz");
                } else if (dataRateCode[0]==(byte)0x03) {
                    DATA_RATE_SAMPLES_PER_SECOND = 1000;
                    Log.e("Data Rate:","1000Hz");
                } else {
                    Log.e("Data Rate:"," ERROR: NOT RECOGNIZED!");
                }
                runOnUiThread(new Runnable() {
                    @Override
                    public void run() {
                        mDataRate.setText("Data Rate: "+String.valueOf(dataRateCode[0]));
                    }
                });*/
            }
        } else {
            Log.e(TAG, "onCharacteristic Read Error" + status);
        }
    }
    private int dataCnt1000 = 0;
    @Override
    public void onCharacteristicChanged(BluetoothGatt gatt, BluetoothGattCharacteristic characteristic) {
        if (AppConstant.CHAR_ION_NA_SIGNAL.equals(characteristic.getUuid())) {
            int dataIonSensor = characteristic.getIntValue(BluetoothGattCharacteristic.FORMAT_UINT8, 0);
//            updateIonSensorState(dataIonSensor);
        }
        //TODO: ADD BATTERY MEASURE CAPABILITY IN FIRMWARE: (ble_ADC)
        if (AppConstant.CHAR_BATTERY_LEVEL.equals(characteristic.getUuid())) {
            batteryLevel = characteristic.getIntValue(BluetoothGattCharacteristic.FORMAT_UINT8, 0);
            updateBatteryStatus(batteryLevel, batteryLevel + " %");
            String timeStamp = getTimeStamp2();
            exportLogFile(false, "Battery Level Changed at " + timeStamp + getDetails()+"\r\n");
            /*
            try {
                GmailSender sender = new GmailSender("developmenttestingmsm@gmail.com",">S#AWr+L?ZwBQ'y");
                sender.sendMail("ECG Device - Battery Level Update","Level Changed at "+timeStamp+getDetails(),"developmenttestingmsm@gmail.com","musasmahmood@gmail.com");
            } catch (Exception e) {
                Log.e("SendMail", e.getMessage(), e);
            }*/
            Log.i(TAG, "Battery Level :: " + batteryLevel);
        }
        if (AppConstant.CHAR_EMG_SIGNAL.equals(characteristic.getUuid())) {
            byte[] dataEmgBytes = characteristic.getValue();
            int byteLength = dataEmgBytes.length;
            //TODO: Remember to check/uncheck plotImplicitXVals (boolean)
            if(!plotImplicitXVals) {
                if(DATA_RATE_SAMPLES_PER_SECOND==250) {
                    for (int i = 0; i < byteLength/2; i++) { //0→9
                        dataCnt1000++;
                        writeToDisk(( unsignedToSigned(unsignedBytesToInt(dataEmgBytes[2*i],dataEmgBytes[2*i+1]),16) ));
                        shortUpdateECGState( unsignedToSigned(unsignedBytesToInt(dataEmgBytes[2*i],dataEmgBytes[2*i+1]),16) );
                    }
                } else if (DATA_RATE_SAMPLES_PER_SECOND==500) {
                    for (int i = 0; i < byteLength/2; i++) { //0→9
                        writeToDisk(( unsignedToSigned(unsignedBytesToInt(dataEmgBytes[2*i],dataEmgBytes[2*i+1]),16) ));
                        dataCnt1000++;
                        if((i+1)%2==0) {
                            shortUpdateECGState( unsignedToSigned(unsignedBytesToInt(dataEmgBytes[2*i],dataEmgBytes[2*i+1]),16) );
                        }
                    }
                    // Add this later, will fir_combined_500 + own graphing method
                } else if (DATA_RATE_SAMPLES_PER_SECOND==1000) {
                    for (int i = 0; i < byteLength/2; i++) {
                        double c = ( unsignedToSigned(unsignedBytesToInt(dataEmgBytes[2*i],dataEmgBytes[2*i+1]),16) );
//                        writeToDisk(( unsignedToSigned(unsignedBytesToInt(dataEmgBytes[2*i],dataEmgBytes[2*i+1]),16) ));
                        writeToDisk((int)c);
                        final double cc = 2.42*c/32767.0;
                        dataCnt1000++;
                        if(dataCnt1000%4==0) {
                            runOnUiThread(new Runnable() {
                                @Override
                                public void run() {
                                    mEcgValues.setText(" " + " Voltage="+String.format("%1.4f",cc)+"V");
                                }
                            });

                            shortUpdateECGState( unsignedToSigned(unsignedBytesToInt(dataEmgBytes[2*i],dataEmgBytes[2*i+1]),16) );
                        } else if (dataCnt1000 == 2147483644+1) {// max value
                            dataCnt1000=0;
                        }
                    }
                }
            }
//            Log.e("Data = ",String.valueOf(dataEmgBytes.length)+" # of bytes");
            getDataRate(dataEmgBytes.length / 2);
        }
        if (AppConstant.CHAR_MPU_COMBINED.equals(characteristic.getUuid())) {
            byte[] allByteData = characteristic.getValue();
            mpuData = getIntArrFromByteArr(allByteData);
            //LOADS ACC/GYRO DATA TO END OF ARRAY SECTION, e.g. data[99]=AccX; data[199]=AccY
            //...data[599]=GyroZ.
            for (int i = 0; i < 6; i++) {
                allMpuData[(i + 1) * 100 - 1] = (double) mpuData[i];
            }
            //SHIFT EVERYTHING BACKWARDS BY 1 (LAST VAL IN EACH SET REMAINS SAME, e.g. 98&99):
            for (int i = 0; i < (100 - 1); i++) {
                for (int j = 0; j < 6; j++) {
                    allMpuData[i + j * 100] = allMpuData[1 + (i + j * 100)];
                }
            }
            double[] returned = new double[8];
            //jTestMasterNew3(allMpuData); //should be size 8;
            //TODO: Replace All Manual Array Copies!
            for (int i = 0; i < 5; i++) {
                yfit[i] = returned[i];
            }
            //TODO: e.g System.arraycopy(returned, 0, yfit, 0, yfit.length);
            accelerometer_resultant = returned[5];
//            updateReturnedMpuData(yfit, accelerometer_resultant);
            getDataRate2();
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
    private double[] ECGBufferAfibDetection = new double[3000];

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
                if(ecgDataSeries.size()>249) {
                    double max = findGraphMax(ecgDataSeries);
                    double min = findGraphMin(ecgDataSeries);
                    if((max-min)!=0) {
                        if(currentBM!=BoundaryMode.AUTO) {
                            ecgPlot.setRangeBoundaries(-2.5, 2.5, BoundaryMode.AUTO);
                            currentBM = BoundaryMode.AUTO;
                        }
                        ecgPlot.setRangeStepValue((max-min)/5);
                    } else {
                        if(currentBM!=BoundaryMode.FIXED) {
                            ecgPlot.setRangeBoundaries(min-1, max+1, BoundaryMode.FIXED);
                            currentBM = BoundaryMode.FIXED;
                        }
                        ecgPlot.setRangeStepValue(2.0/5.0);
                    }
                }
                //TODO: Now plot:
                 //The idea is that every time this thing is triggered, we want to plot the entirety
                 //of [explicitXValsLong, ECGBufferFiltered2] BEFORE the next second occurs.
                if(mTimeSeconds == HISTORY_SECONDS_2) {
                    newMinX = Math.floor(explicitXValsLong[0]);
                    newMaxX = Math.floor(explicitXValsLong[999]);
                    ecgPlot.setDomainBoundaries(newMinX, newMaxX, BoundaryMode.AUTO);
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
            mEcgValues.setText(" " + " Voltage="+String.format("%1.4f",dataVoltage)+"V");
            }
        });
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

    class graphEcgData3 implements Runnable {

        @Override
        public void run() {
            if (startedGraphingEcgData) {
                if (ecgIndex < 1000) {
                    if (ecgDataSeries.size() > 1000) {
                        ecgDataSeries.removeFirst();
                    }
                    if (plotImplicitXVals) {
                        if (filterData) {
                            ecgDataSeries.addLast(null, ECGBufferFiltered2[ecgIndex]);
                            /*if(!filterType){
                            } else {
                                ecgDataSeries.addLast(null, ECGBufferFiltered3[ecgIndex]);
                            }*/
                        } else {
                            ecgDataSeries.addLast(null, ECGBufferUnfiltered[ecgIndex]);
                        }
                    } else {
                        if (filterData) {
                            ecgDataSeries.addLast(explicitXValsLong[ecgIndex], ECGBufferFiltered2[ecgIndex]);
                            /*if(!filterType){
                            } else {
                                ecgDataSeries.addLast(explicitXValsLong[ecgIndex], ECGBufferFiltered3[ecgIndex]);
                            }*/
                        } else {
                            ecgDataSeries.addLast(explicitXValsLong[ecgIndex], ECGBufferUnfiltered[ecgIndex]);
                        }
                    }
                    ecgIndex++;
                } else {
                    ecgIndex = 750;
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

    private static final char[] HEX_CHARS = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A',
            'B', 'C', 'D', 'E', 'F' };
    public static String toHexString(byte[] bytes) {
        char[] hexChars = new char[bytes.length * 2];
        int v;
        for (int j = 0; j < bytes.length; j++) {
            v = bytes[j] & 0xFF;
            hexChars[j * 2] = HEX_CHARS[v >>> 4];
            hexChars[j * 2 + 1] = HEX_CHARS[v & 0x0F];
//            hexChars[j * 2 + 1] = HEX_CHARS[v >>> 4];
//            hexChars[j * 2] = HEX_CHARS[v & 0x0F];
        }
        return new String(hexChars);
    }

//    String toBinary( byte[] bytes ) {
//        StringBuilder sb = new StringBuilder(bytes.length * Byte.SIZE);
//        for( int i = 0; i < Byte.SIZE * bytes.length; i++ ){
//            sb.append((bytes[i / Byte.SIZE] << i % Byte.SIZE & 0x80) == 0 ? '0' : '1');
//            if((i+1)%8==0 && i>0)sb.append(" ");
//        }
//        return sb.toString();
//    }

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
//                        mDataRate.setText("0 Hz");
                    }
                });
                exportLogFile(false, "Connected @ "+getTimeStamp()+"\r\n"+"\r\n");

                /*String timeStamp = getTimeStamp2();
                try {
                    GmailSender sender = new GmailSender("developmenttestingmsm@gmail.com",">S#AWr+L?ZwBQ'y");
//                    GmailSender sender = new GmailSender("halamadridgoals@gmail.com","Axjy1z77");
                    sender.sendMail("ECG Device Connected","Connected at "+timeStamp+getDetails(),"developmenttestingmsm@gmail.com","musasmahmood@gmail.com");
//                    sender.sendMail("ECG Device Connected","Connected at "+timeStamp+getDetails(),"halamadridgoals@gmail.com","musasmahmood@gmail.com");
                } catch (Exception e) {
                    Log.e("SendMail", e.getMessage(), e);
                }*/
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

                /*String timeStamp2 = getTimeStamp2();
                try {
                    GmailSender sender = new GmailSender("developmenttestingmsm@gmail.com",">S#AWr+L?ZwBQ'y");
//                    GmailSender sender = new GmailSender("halamadridgoals@gmail.com","Axjy1z77");
                    sender.sendMail("ECG Device Disconnected","Disconnected at "+timeStamp2+getDetails(),"developmenttestingmsm@gmail.com","musasmahmood@gmail.com");
//                    sender.sendMail("ECG Device Disconnected","Disconnected at "+timeStamp2+getDetails(),"halamadridgoals@gmail.com","musasmahmood@gmail.com");
                } catch (Exception e) {
                    Log.e("SendMail", e.getMessage(), e);
                }*/
                //TODO: ATTEMPT TO RECONNECT:
                /*if (mBluetoothLe != null) {
                    mBluetoothLe.connect(mBluetoothDevice, false);
                } connect();*/
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
                if (mBluetoothGatt == null /*||
                        mBluetoothAdapter == null*/ ||
                        !mConnected) {
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
//                mBatteryLevel.setText("Not Available");
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
    private native int jmainFirFilter(boolean b);

    private native double[] jFirFilter(double[] ecg); //size = 1000

    //ECG BW Filter:

    private native double[] jBwFilter(double[] ecg);

    //ECG ANALYSIS (SIMPLE):
//    private native double jEcgAnalysis(double[] ecg);

    //ECG Afib Detection:
//    private native int jmainAfibDetectionInit(boolean b);

//    private native double jAfibDetectionAnalysis(double[] ecg_h, double MeanHeartRate);

    //FALL DETECTION:
}
