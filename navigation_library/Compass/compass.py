import numpy as np
from scipy.optimize import curve_fit


def fit_f(x, p1, p2, p3):
    # function for approximation of X-gyro output modulated by Earth rotation rate
    # OmegaX = Omega0 + ScaleFactor * EarthRate * cos(yaw + initial_azimuth)
    # here: p1 = Omega0 - gyro bias
    #       p2 = ScaleFactor - gyro scale factor
    #       p3 - initial azimuth
    #       x - current yaw angle
    return p1 + p2 * np.cos((p3 + x) * np.pi / 180)


class Caruseling:
    # class for azimuth estimation using method of X-gyro modulation with Earth rotation rate and reconstruction of
    # initial azimuth angle
    # Each time gyro should be oriented at new yaw angle, stay still in this position (with measuring and smoothing
    # X-gyro output for about 1-2min. Then feed new pair od data (yaw, gyro_out) to class Caruseling and get
    # new initial azimuth estimation. Estimation is provided using the Least Square Method.
    def __init__(self):
        self.yaw = []
        self.gyro_out = []
        self.azimuth = None

    def __call__(self, yaw_in_nav_frame, gyro_out):
        # inputs - next pair of yaw and X-gyro measurement
        # output - estimated initial azimuth of Navigation Frame
        self.yaw.append(yaw_in_nav_frame)
        self.gyro_out.append(gyro_out)
        if len(self.yaw) > 3:
            x = np.array(self.yaw)
            y = np.array(self.gyro_out)
            popt = curve_fit(fit_f, x, y, p0=(0.0, 10.2, 0))
            azimuth = popt[2]
            if popt[1] < 0:
                azimuth += 180
        return self.azimuth


class Maytagging:
    # class for azimuth estimation using maytagging method:
    # fast estimation can be done orienting X-gyro 4 perpendicular directions (four_point method)
    # at each point (direction) gyro should be stable and its output should be smoothed during 1-2min
    # Measurements in 4-points method are:
    # Omega0 = Omega0 + ScaleFactor * EarthRate * cos(azimuth)
    # Omega90 = Omega0 - ScaleFactor * EarthRate * sin(azimuth)
    # Omega180 = Omega0 - ScaleFactor * EarthRate * cos(azimuth)
    # Omega270 = Omega0 + ScaleFactor * EarthRate * sin(azimuth)
    # so:
    #   Omega0 - Omega180 = 2 * EarthRate * cos(azimuth)
    #   Omega270 - Omega90 = 2 * EarthRate * sin(azimuth)
    # and: tg(azimuth) = (Omega270 - Omega90) / (Omega0 - Omega180)
    # two_point method is less accurate due to lack of precise information of EarthRate in current location
    # It uses the fact that Omega0 - Omega180 = 2 * EarthRate * sin(azimuth)
    # for best SNR Omega0 and Omega180 should be measured in approx East (Omega270 in 4-points method)
    # and West (Omega90 in 4-points method) directions.

    def __init__(self):
        self.azimuth = 0

    def two_point(self, gyro_out_at_yaw_0deg, gyro_out_at_yaw_180deg, earth_rate_mod):
        # two X-gyro measurements - in aprrox East direction (0deg) and approx West direction (180deg)
        # output - estimated initial azimuth of Navigation Frame
        azimuth = np.arsin((gyro_out_at_yaw_0deg - gyro_out_at_yaw_180deg) / 2 / earth_rate_mod)
        self.azimuth = (self.azimuth + azimuth) / 2

    def four_point(self, gyro_out_at_yaw_0_90_180_270):
        # input - list of 4 measurements of X-gyro output at 4 perpendicular directions
        # output - estimated initial azimuth of Navigation Frame
        azimuth = np.arctan2(-gyro_out_at_yaw_0_90_180_270[1] - gyro_out_at_yaw_0_90_180_270[3],
                             gyro_out_at_yaw_0_90_180_270[0] - gyro_out_at_yaw_0_90_180_270[2])
        self.azimuth = (self.azimuth + azimuth) / 2
