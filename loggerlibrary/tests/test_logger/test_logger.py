import time
from unittest import TestCase

import h5py

from loggerlibrary.logger.logger import Logger, LoggerType

text_results = ''
bin_results = ''


class TestLogger(TestCase):
    def test_text_logger(self):
        self.logger = Logger('test')
        self.logger.start()

        for i in range(11):
            self.logger('test' + str(i))

        self.logger.stop()
        self.logger.join()

        with open('Logs/test.txt') as fid:
            data = fid.readlines()

        self.assertEqual(len(data), 11)
        self.assertEqual(''.join(data).count('test'), 11)

    def test_text_logger_performance(self):
        self.logger = Logger('test')
        self.logger.start()
        N = 1000000
        txt = 't' * 100

        t0 = time.time()
        for i in range(N):
            self.logger(txt)

        self.logger.stop()
        self.logger.join()
        t1 = time.time()
        global text_results
        text_results = '\n\nLogger text performance: {:0.3f}s for 1M records of 100 chars,' \
                       ' {:0.1f}k records/s\n\n'.format(t1 - t0, 1e6 / (t1 - t0) / 1e3)

    def test_binary_logger(self):
        self.logger = Logger('test_bin', logger_type=LoggerType.Yaw)
        self.logger.start()
        err = False
        N = 2001

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

        with h5py.File('Logs/test_bin.df', 'r') as fid:
            keys = fid.keys()
            self.assertEqual('t' in keys, True)
            self.assertEqual('yaw' in keys, True)
            self.assertEqual('goal_yaw' in keys, True)
            data = fid['yaw']
            self.assertEqual(list(data), [2 * i for i in range(N)])

    def test_binary_logger_performance(self):
        self.logger = Logger('test_bin', logger_type=LoggerType.Velocity)
        self.logger.start()
        err = False
        N = 1000000

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
        global bin_results
        bin_results = '\n\nLogger binary (H5DF) performance: {:0.3f}s for 1M records of 5 fields, ' \
                      '{:0.1f}k records/s\n\n'.format(t1 - t0, 1e6 / (t1 - t0) / 1e3)

    def tearDown(self):
        global text_results
        global bin_results
        if len(text_results) > 0:
            print(text_results)
            text_results = ''
        if len(bin_results) > 0:
            print(bin_results)
            bin_results = ''
