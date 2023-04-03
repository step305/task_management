import requests
import random
import string
import json
from SC_encryption import SC_encryption

class big_brother():
    def __init__(self, address, keys_file="keys.json", series="NumLock"):
        self.address = address
        self.keys_file = keys_file
        try:
            with open(keys_file) as f:
                d = json.load(f)
            self.robot = d['robot']
            self.series = series
            self.encryption = SC_encryption()
            self.encryption.set_jwt_key(d["jwt_key"])
            self.encryption.set_key(self.encryption.decode_key(d["key"]))
        except FileNotFoundError:
            self.encryption = SC_encryption()
            self.register(keys_file, series)
            self.series = series
        return

    def register(self, keys_file, series):
        r = requests.post(self.address+'get_sn/', data={
            'series': series,
        })
        try:
            self.robot = r.json()['SN']
        except:
            print(r.json()['Result'])
            exit()
        new_key = r.json()['key']
        new_jwt_key = r.json()['jwt_key']
        with open(self.keys_file, 'w+') as f:
            d = {
                "robot": self.robot,
                "series": series,
                "key": new_key,
                "jwt_key": new_jwt_key
            }
            json.dump(d, f)
        self.encryption.set_jwt_key(new_jwt_key)
        self.encryption.set_key(self.encryption.decode_key(new_key))
        return r

    def new_session(self):
        encoded = self.encryption.encode_jwt({'data':self.robot})
        r = requests.post(self.address+'start_session/', data={
            'data': encoded,
        })
        new_key = self.encryption.decode_data(r.json()['data'])['key']
        self.encryption.set_key(self.encryption.decode_key(new_key))
        with open(self.keys_file, 'w+') as f:
            d = {
                "robot": self.robot,
                "series": self.series,
                "key": new_key,
                "jwt_key": self.encryption.get_jwt_key()
            }
            json.dump(d, f)
        return r

    def report_sensor(self, sensor, node, topic, data):
        encoded = self.encryption.code_data({
            'sensor': sensor,
            'node': node,
            'topic': topic,
            'data': data,
        })
        encoded = self.encryption.encode_jwt({'data':encoded})
        r = requests.post(self.address+'report_sensor/', data={
            'robot': self.robot,
            'data': encoded,
        })
        return r

    def session_file_upload(self, bag_file_name=False, media_file_name=False, GPS_file_name=False):
        files = []
        encoded = self.encryption.encode_jwt({'data':self.robot})
        if bag_file_name:
            files.append(('bag', open(bag_file_name, 'rb')))
        if media_file_name:
            files.append(('media', open(media_file_name, 'rb')))
        if GPS_file_name:
            files.append(('GPS', open(GPS_file_name, 'rb')))
        r = requests.post(self.address+'session_upload/', data={'data': encoded}, files=files)
        return r

    def report_session(self, av_lon=0, av_lat=0, volt_start=-1, volt_end=-1, GPS_status=-1):
        encoded = self.encryption.code_data({
            'av_lon': av_lon,
            'av_lat': av_lat,
            'volt_start': volt_start,
            'volt_end': volt_end,
            'GPS_status': GPS_status
        })
        encoded = self.encryption.encode_jwt({'data':encoded})
        r = requests.post(self.address+'report_session/', data={
            'robot': self.robot,
            'data': encoded,
        })
        return r

    def close_session(self):
        encoded = self.encryption.encode_jwt({'data':self.robot})
        r = requests.post(self.address+'close_session/', data={
            'data': encoded,
        })
        return r

    def get_firmware_version(self, type):
        encoded = self.encryption.encode_jwt({'data':self.robot, 'type': type})
        r = requests.post(self.address+'get_latest_firmware_version/', data={
            'data': encoded,
        })
        result = r.json()
        return result['Result']


    def get_firmware(self, type, version, fname):
        encoded = self.encryption.encode_jwt({'data':self.robot, 'type': type, 'version': version})
        r = requests.post(self.address+'get_latest_firmware/', data={
            'data': encoded
        })
        handle = open(fname, "wb")
        for chunk in r.iter_content(chunk_size=512):
            if chunk:  # filter out keep-alive new chunks
                handle.write(chunk)
        handle.close()
        return r
