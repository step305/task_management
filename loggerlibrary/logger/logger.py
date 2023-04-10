import os
import threading
import queue
import datetime
import time
import enum
import h5py

BUF_LEN_BINARY = 1000
BUF_LEN_TEXT = 10


class LoggerType(enum.Enum):
    Text = ''
    Coords = ('t', 'x', 'y', 'goal_x', 'goal_y')
    Velocity = ('t', 'vx', 'vy', 'goal_vx', 'goal_vy',)
    Yaw = ('t', 'yaw', 'goal_yaw')
    Rate = ('t', 'wx', 'wy', 'wz', 'goal_wz')
    Accel = ('t', 'ax', 'ay', 'az')
    Odometer = ('t', 'l_odo', 'r_odo', 'l_goal_odo', 'r_goal_odo')
    Magnetometer = ('t', 'mx', 'my', 'mz')
    Control = ('t', 'l_control_wheel', 'r_control_wheel')
    Battery = ('t', 'voltage', 'charge')


class Logger(threading.Thread):
    def __init__(self, logger_name, logger_type=LoggerType.Text):
        super(Logger, self).__init__()
        self.stop_event = threading.Event()
        self.stop_event.clear()
        self.input_queue = queue.Queue(1000)
        self.busy = False
        if not os.path.exists('Logs'):
            os.mkdir('Logs')
        self.out_file_name = os.path.join('Logs', logger_name)
        self.start_log = datetime.datetime.now()
        self.logger_type = logger_type
        if self.logger_type == LoggerType.Text:
            self.file_id = open(self.out_file_name + '.txt', 'w')

        else:
            self.file_id = h5py.File(self.out_file_name + '.df', 'w')

    def run(self):
        if not self.logger_type == LoggerType.Text:
            N = len(self.logger_type.value)
            for v in self.logger_type.value:
                self.file_id.create_dataset(v,
                                            (0,),
                                            dtype='<f4',
                                            chunks=(1000,),
                                            compression='gzip',
                                            maxshape=(None,)
                                            )
            buffer = {v: [] for v in self.logger_type.value}
        else:
            buffer = []
        while not self.stop_event.is_set():
            self.busy = True
            dt, data = self.input_queue.get(block=True)
            if dt is None:
                break
            self.busy = False
            if self.logger_type == LoggerType.Text:
                buffer.append('{:02d}:{:02d}:{:02d}.{:d} {:02d}.{:02d}.{:04d}: {}'.format(
                    dt.hour, dt.minute, dt.second, dt.microsecond, dt.day, dt.month, dt.year, data
                ))
                if len(buffer) == BUF_LEN_TEXT:
                    self.file_id.write('\n'.join(buffer) + '\n')
                    self.file_id.flush()
                    buffer = []
            else:
                for i, v in enumerate(self.logger_type.value):
                    buffer[v].append(data[i])
                if len(buffer[self.logger_type.value[0]]) >= BUF_LEN_BINARY:
                    for v in self.logger_type.value:
                        self.file_id[v].resize((self.file_id[v].shape[0] + len(buffer[v]),))
                        self.file_id[v][-len(buffer[v]):] = buffer[v]
                    buffer = {v: [] for v in self.logger_type.value}
                    self.file_id.flush()
        while not self.input_queue.empty():
            dt, data = self.input_queue.get()
            if dt is None:
                break
            if self.logger_type == LoggerType.Text:
                buffer.append('{:02d}:{:02d}:{:02d}.{:d} {:02d}.{:02d}.{:04d}: {}'.format(
                    dt.hour, dt.minute, dt.second, dt.microsecond, dt.day, dt.month, dt.year, data
                ))
            else:
                for i, v in enumerate(self.logger_type.value):
                    buffer[v].append(data[i])
        if self.logger_type == LoggerType.Text:
            if len(buffer) > 0:
                self.file_id.write('\n'.join(buffer) + '\n')
        else:
            if len(buffer[self.logger_type.value[0]]) > 0:
                for v in self.logger_type.value:
                    self.file_id[v].resize((self.file_id[v].shape[0] + len(buffer[v]),))
                    self.file_id[v][-len(buffer[v]):] = buffer[v]
        self.file_id.flush()
        self.file_id.close()

    def __call__(self, data):
        if not self.logger_type == LoggerType.Text:
            if not len(data) == len(self.logger_type.value):
                self.stop()
                raise AssertionError('Wrong number of values in logger {}: {} given but {} expected'.
                                     format(self.out_file_name,
                                            len(data),
                                            len(self.logger_type.value)))
        self.input_queue.put((datetime.datetime.now(), data))

    def stop(self):
        self.stop_event.set()
        if self.busy:
            self.input_queue.put((None, None))
        time.sleep(0.1)
