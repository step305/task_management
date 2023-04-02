# Compass
## Description

***Compass*** is a library for azimuth estimation using two methods - caruseling or maytagging.<br>
It depends on externally provided data on azimuth gyro output (in deg/hr) and yaw in navigation frame.<br>
Output is estimated azimuth relative to navigation frame.

## Usage

Before use instance of ***Caruseling*** or ***Maytagging*** should be created:
```python
from compass import Caruseling, Maytagging
# for caruseling method
estimator_caruseling = Caruseling()
# for maytagging method
estimator_maytagging = Maytagging()

# estimation
yaw = 15  # in degrees
gyro_out = 10  # in deg/hr
# get gyro out and current yaw somewhere
azimuth = estimator_caruseling(yaw, gyro_out)
# if not ready azimuth is equal to None

# get gyro_out [deg/hr] at yaw = 0, 90, 180, 270 degrees
gyro_out = [1, 11, 1, -9]
earth_rate = 10  # deg/hr 
azimuth = estimator_maytagging.two_point(gyro_out[0], gyro_out[2], earth_rate)
# or
azimuth = estimator_maytagging.four_point(gyro_out)
```

## Testing

Omitted due to trivial nature.
