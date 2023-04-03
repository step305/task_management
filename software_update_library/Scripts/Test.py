
from big_brother import big_brother


s = big_brother('http://127.0.0.1:8000/RobotMonitor/', 'keys.json')

version = s.get_firmware_version('ROS')
r = s.get_firmware('ROS', '1.0.0', 'archieve.run')

# r = s.new_session()
# r = s.session_file_upload(media_file_name="images.png")
# r = s.close_session()

print(version)
# print(r)
