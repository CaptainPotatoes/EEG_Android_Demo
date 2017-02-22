package com.mahmoodms.bluetooth.ecgfallsensordemo;

import android.Manifest;
import android.app.ActionBar;
import android.app.Activity;
import android.bluetooth.BluetoothAdapter;
import android.bluetooth.BluetoothManager;
import android.bluetooth.le.ScanCallback;
import android.bluetooth.le.ScanResult;
import android.content.Context;
import android.content.Intent;
import android.content.pm.ActivityInfo;
import android.content.pm.PackageManager;
import android.graphics.Color;
import android.graphics.drawable.ColorDrawable;
import android.os.Bundle;
import android.os.Handler;
import android.support.v4.app.ActivityCompat;
import android.support.v4.content.ContextCompat;
import android.view.Menu;
import android.view.MenuItem;
import android.view.View;
import android.view.WindowManager;
import android.widget.AdapterView;
import android.widget.ListView;
import android.widget.Toast;

import java.util.ArrayList;

/**
 * Created by mahmoodms on 6/30/2016.
 */



public class MainActivity extends Activity {
    private boolean mScanning;
    private Handler mHandler;
    private ListView scanningDeviceListView;
    private BluetoothAdapter mBluetoothAdapter;
    private static final int MY_PERMISSIONS_LOCATIONS_COARSE = 156;
    private static final int MY_PERMISSIONS_WRITE_STORAGE = 754;
    private static final int MY_PERMISSIONS_SEND_SMS = 2285;
    private static final int MY_PERMISSIONS_INTERNET = 1654;
    private static final int MULTIPLE_PERMISSIONS_REQUEST = 139;

    private ScannedDeviceAdapter scannedDeviceAdapter;
    private static final long SCAN_PERIOD = 10000;
    private final static int REQUEST_ENABLE_BT = 12;

    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_device_scan);
        //Set Orientation in Landscape
        setRequestedOrientation(ActivityInfo.SCREEN_ORIENTATION_LANDSCAPE);
        //Set up action bar:
        if (getActionBar()!=null) getActionBar().setDisplayHomeAsUpEnabled(false);
        ActionBar actionBar = getActionBar();
        actionBar.setBackgroundDrawable(new ColorDrawable(Color.parseColor("#800020")));
        //Set orientation of device based on screen type/size:
        //Flag to keep screen on (stay-awake):
        getWindow().addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);
        //Set up timer handler:
        mHandler = new Handler();

        //Initialize scanningDeviceListView Adapter:
        scanningDeviceListView = (ListView) findViewById(R.id.scanningList);

        //Check for BLE Support
        if (!getPackageManager().hasSystemFeature(PackageManager.FEATURE_BLUETOOTH_LE)) {
            Toast.makeText(this, R.string.ble_not_supported, Toast.LENGTH_SHORT).show();
            finish();
        }
        //Initialize Bluetooth manager
        BluetoothManager bluetoothManager = (BluetoothManager) getSystemService(Context.BLUETOOTH_SERVICE);
        mBluetoothAdapter = bluetoothManager.getAdapter();
        if(mBluetoothAdapter == null) {
            Toast.makeText(this,R.string.error_bluetooth_not_supported,Toast.LENGTH_SHORT).show();
            finish();
            return;
        }
        /**
         * TODO: Need to add request that will give option to enable permissions, if permissions have not been granted yet
         */
        int permissionCheck = ContextCompat.checkSelfPermission(MainActivity.this, Manifest.permission.ACCESS_COARSE_LOCATION);
        int permissionCheck2 = ContextCompat.checkSelfPermission(MainActivity.this, Manifest.permission.WRITE_EXTERNAL_STORAGE);
        int permissionCheck3 = ContextCompat.checkSelfPermission(MainActivity.this, Manifest.permission.SEND_SMS);
        int permissionCheck4 = ContextCompat.checkSelfPermission(MainActivity.this, Manifest.permission.INTERNET);
        if(permissionCheck!=PackageManager.PERMISSION_GRANTED && permissionCheck2!=PackageManager.PERMISSION_GRANTED &&  permissionCheck3!=PackageManager.PERMISSION_GRANTED) {
            ActivityCompat.requestPermissions(MainActivity.this, new String[]{Manifest.permission.ACCESS_COARSE_LOCATION,
                    Manifest.permission.WRITE_EXTERNAL_STORAGE, Manifest.permission.SEND_SMS, Manifest.permission.INTERNET}, MULTIPLE_PERMISSIONS_REQUEST);
        } else if (permissionCheck!=PackageManager.PERMISSION_GRANTED) {
            ActivityCompat.requestPermissions(MainActivity.this, new String[]{Manifest.permission.ACCESS_COARSE_LOCATION}, MY_PERMISSIONS_LOCATIONS_COARSE);
        } else if (permissionCheck2!=PackageManager.PERMISSION_GRANTED) {
            ActivityCompat.requestPermissions(MainActivity.this, new String[]{Manifest.permission.WRITE_EXTERNAL_STORAGE}, MY_PERMISSIONS_WRITE_STORAGE);
        } else if (permissionCheck3!=PackageManager.PERMISSION_GRANTED) {
            ActivityCompat.requestPermissions(MainActivity.this, new String[]{Manifest.permission.SEND_SMS}, MY_PERMISSIONS_SEND_SMS);
        } else if (permissionCheck4!=PackageManager.PERMISSION_GRANTED) {
            ActivityCompat.requestPermissions(MainActivity.this, new String[]{Manifest.permission.INTERNET}, MY_PERMISSIONS_INTERNET);
        }

        //Initialize list view adapter:
        scannedDeviceAdapter = new ScannedDeviceAdapter(this, R.layout.scanning_item, new ArrayList<ScannedDevice>());
        scanningDeviceListView.setAdapter(scannedDeviceAdapter);
        // Click Item Listener:
        scanningDeviceListView.setOnItemClickListener(new AdapterView.OnItemClickListener() {
            @Override
            public void onItemClick(AdapterView<?> adapterView, View view, int position, long id) {
            ScannedDevice item = scannedDeviceAdapter.getItem(position);
            if (item!=null) {
                final Intent intent = new Intent(MainActivity.this, DeviceControlActivity.class);
                intent.putExtra(AppConstant.EXTRAS_DEVICE_NAME, item.getDisplayName());
                intent.putExtra(AppConstant.EXTRAS_DEVICE_ADDRESS, item.getDeviceMac());
                if(mScanning) {
                    //TODO - Requires v21 or greater...
                    //@Deprecated:
//                    mBluetoothAdapter.stopLeScan(mLeScanCallback);
                    mBluetoothAdapter.getBluetoothLeScanner().stopScan(mScanCallback);
                    mScanning = false;
                }
                startActivity(intent);
            }
            }
        });
    }

    @Override
    public void onRequestPermissionsResult(int requestCode, String permissions[], int[] grantResults) {
        switch (requestCode) {
            case MY_PERMISSIONS_LOCATIONS_COARSE: {
                if(grantResults.length>0 && grantResults[0]==PackageManager.PERMISSION_GRANTED) {

                } else {
                    runOnUiThread(new Runnable() {
                        @Override
                        public void run() {
//                            Toast toast = Toast.makeText(getApplicationContext(),
//                                    "This permission is needed for Bluetooth LE!", Toast.LENGTH_LONG);
//                            toast.show();
                        }
                    });
                }
                return;
            }
            case MY_PERMISSIONS_WRITE_STORAGE: {
                if(grantResults.length>0 && grantResults[0]==PackageManager.PERMISSION_GRANTED) {

                } else {
                    runOnUiThread(new Runnable() {
                        @Override
                        public void run() {
//                            Toast toast = Toast.makeText(getApplicationContext(),
//                                    "This permission is required to save data!", Toast.LENGTH_LONG);
//                            toast.show();
                        }
                    });
                }
                return;
            }
            case MULTIPLE_PERMISSIONS_REQUEST: {
                if(grantResults.length>0 && grantResults[0]==PackageManager.PERMISSION_GRANTED && grantResults[1]==PackageManager.PERMISSION_GRANTED) {
                } else {
                    runOnUiThread(new Runnable() {
                        @Override
                        public void run() {
                            Toast toast = Toast.makeText(getApplicationContext(),
                                    "These permissions are required!", Toast.LENGTH_LONG);
                            toast.show();
                        }
                    });
                }
            }
        }
    }

    @Override
    protected void onPause() {
        super.onPause();
        scanLeDevice(false);
        scannedDeviceAdapter.clear();
    }

    @Override
    protected void onResume() {
        super.onResume();
        /**
         * Ensures Bluetooth is enabled on the device - if not enabled - fire intent to display a
         * dialog to ask permission to enable
         */
        if (!mBluetoothAdapter.isEnabled()) {
            if (!mBluetoothAdapter.isEnabled()) {
                Intent enableBt = new Intent(BluetoothAdapter.ACTION_REQUEST_ENABLE);
                startActivityForResult(enableBt, REQUEST_ENABLE_BT);
            }
        }
        scanLeDevice(true);
    }

    private ScanCallback mScanCallback = new ScanCallback() {
        @Override
        public void onScanResult(int callbackType, final ScanResult result) {
            runOnUiThread(new Runnable() {
                @Override
                public void run() {
                    scannedDeviceAdapter.update(result.getDevice(), result.getRssi(), result.getScanRecord());
                    scannedDeviceAdapter.notifyDataSetChanged();
                }
            });
            super.onScanResult(callbackType, result);
        }
    };

    /*private BluetoothAdapter.LeScanCallback mLeScanCallback = new BluetoothAdapter.LeScanCallback() {
        @Override
        public void onLeScan(final BluetoothDevice device,
                             final int rssi,
                             final byte[] scanRecord) {
            runOnUiThread(new Runnable() {
                @Override
                public void run() {
                    scannedDeviceAdapter.update(device, rssi, scanRecord);
                    scannedDeviceAdapter.notifyDataSetChanged();
                }
            });
        }
    };
*/
    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
        getMenuInflater().inflate(R.menu.menu_device_scan, menu);
        if (!mScanning) {
            menu.findItem(R.id.menu_stop).setVisible(false);
            menu.findItem(R.id.menu_scan).setVisible(true);
            menu.findItem(R.id.menu_refresh).setActionView(null);
        } else {
            menu.findItem(R.id.menu_stop).setVisible(true);
            menu.findItem(R.id.menu_scan).setVisible(false);
            menu.findItem(R.id.menu_refresh).setActionView(
                    R.layout.actionbar_progress);
        }
        return true;
    }

    private void scanLeDevice(final boolean enable) {
        if(enable) {
            //stops scanning after ~seconds
            mHandler.postDelayed(new Runnable() {
                @Override
                public void run() {
                    mScanning = false;
//                    mBluetoothAdapter.stopLeScan(mLeScanCallback);
                    mBluetoothAdapter.getBluetoothLeScanner().stopScan(mScanCallback);
                    invalidateOptionsMenu();
                }
            }, SCAN_PERIOD);
            mScanning = true;
//            mBluetoothAdapter.startLeScan(mLeScanCallback);
            mBluetoothAdapter.getBluetoothLeScanner().startScan(mScanCallback);
        } else {
            mScanning = false;
//            mBluetoothAdapter.stopLeScan(mLeScanCallback);
            mBluetoothAdapter.getBluetoothLeScanner().stopScan(mScanCallback);
        }
        invalidateOptionsMenu();
    }

    @Override
    protected void onActivityResult(int requestCode, int resultCode, Intent data) {
        //if user chose not to enable BT
        if (requestCode == REQUEST_ENABLE_BT && resultCode == Activity.RESULT_CANCELED) {
            finish();
            return;
        }
        super.onActivityResult(requestCode, resultCode, data);
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        switch (item.getItemId()) {
            case R.id.menu_scan:
                scannedDeviceAdapter.clear();
                scanLeDevice(true);
                break;
            case R.id.menu_stop:
                scanLeDevice(false);
                break;
        }
        return true;
    }
}
