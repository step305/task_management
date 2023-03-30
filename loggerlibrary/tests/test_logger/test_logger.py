import time
from unittest import TestCase

import h5py

from loggerlibrary.logger.logger import Logger, LoggerType

text_results = ''
bin_results = ''


class TestLogger(TestCase):
    def test_text_logger(self):
        # test Logger ability to store text data
        self.logger = Logger('test')  # new logger by default is text-logger
        self.logger.start()

        # write some trash in log
        for i in range(11):
            self.logger('test' + str(i))

        # stop and clear logger
        self.logger.stop()
        self.logger.join()

        # read actual data in log-file
        with open('Logs/test.txt') as fid:
            data = fid.readlines()

        # test if actual data is same as written one
        self.assertEqual(len(data), 11)
        self.assertEqual(''.join(data).count('test'), 11)

    def test_text_logger_performance(self):
        # test Logger's performance on text logs
        self.logger = Logger('test')  # new logger by default is text-logger
        self.logger.start()
        N = 1000000  # amount of log data for test
        txt = 't' * 100  # test string to store in log (100 chars)

        # write log and measure time
        t0 = time.time()
        for i in range(N):
            self.logger(txt)

        self.logger.stop()
        self.logger.join()
        t1 = time.time()

        # estimate performance
        global text_results
        text_results = '\n\nLogger text performance: {:0.3f}s for 1M records of 100 chars,' \
                       ' {:0.1f}k records/s\n\n'.format(t1 - t0, 1e6 / (t1 - t0) / 1e3)

    def test_binary_logger(self):
        # test Logger's ability to store binary logs
        # new logger to store yaw measurements (timestamp, yaw, desired yaw)
        self.logger = Logger('test_bin', logger_type=LoggerType.Yaw)
        self.logger.start()
        err = False
        N = 2001  # amount of log data for test, list with 3 numbers will be logged

        # store binary data
        for i in range(N):
            try:
                self.logger([i, 2 * i, i])
            except AssertionError as e:
                print(e)
                err = True
                break

        self.logger.stop()
        try:
            self.logger.join()
        except Exception as e:
            print(e)
        self.assertEqual(err, False)

        # read out and review stored data
        with h5py.File('Logs/test_bin.df', 'r') as fid:  # open H5DF file
            keys = fid.keys()  # read list of keys - names of stored arrays in file
            # test if valid arrays exist: t, yaw, goal_yaw
            self.assertEqual('t' in keys, True)
            self.assertEqual('yaw' in keys, True)
            self.assertEqual('goal_yaw' in keys, True)
            # get data for yaw array
            data = fid['yaw']
            # compare yaw with correct one
            self.assertEqual(list(data), [2 * i for i in range(N)])

    def test_binary_logger_performance(self):
        # test Logger's performance on binary logs
        # new logger to store velocity measurements (timestamp, velx, vely, desired velx, desired vely)
        self.logger = Logger('test_bin', logger_type=LoggerType.Velocity)
        self.logger.start()
        err = False
        N = 1000000  # amount of binary records for test

        # write log and measure time
        # list of 5 values is written
        t0 = time.time()
        for i in range(N):
            try:
                self.logger([i, 2 * i, -i, 3 * i, -2 * i])
            except AssertionError as e:
                print(e)
                err = True
                break

        self.logger.stop()
        try:
            self.logger.join()
        except Exception as e:
            print(e)
        self.assertEqual(err, False)
        t1 = time.time()

        # estimate performance
        global bin_results
        bin_results = '\n\nLogger binary (H5DF) performance: {:0.3f}s for 1M records of 5 fields, ' \
                      '{:0.1f}k records/s\n\n'.format(t1 - t0, 1e6 / (t1 - t0) / 1e3)

    def tearDown(self):
        # utility stuff for nicer output of results
        global text_results
        global bin_results
        if len(text_results) > 0:
            print(text_results)
            text_results = ''
        if len(bin_results) > 0:
            print(bin_results)
            bin_results = ''
