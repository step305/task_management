# Logger
## Description

***Logger*** is a threaded logger.<br>
It can be created for text logs as well as binary logs (in H5DF format).<br>


## Usage

Before use instance of ***Logger*** should be created:
```python
from logger import Logger, LoggerType
# for logging text messages
text_logger = Logger('some-logger-name', logger_type=LoggerType.Text)
# for logging velocity
velocity_logger = Logger('another-name', logger_type=LoggerType.Velocity)

# ... actual logging
text_logger('Some important message')
t = 0.01
vx = vy = 0.3
vx_desired = vy_desired = 4.0
velocity_logger([t, vx, vy, vx_desired, vy_desired])

# ... closing logs
text_logger.stop()
text_logger.join()
velocity_logger.stop()
velocity_logger.join()
```
Logs are created in ***Logs*** sub-directory.<br>
Name of log-file is equal to log's name, extension - '.txt' for text logs and '.df' for binary logs.

## Testing

From main directory (parent of *loggerlibrary*) run:
```shell
python3 -m unittest tests/test_logger/test_logger.py
```

