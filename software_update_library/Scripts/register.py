import requests
import json

class register():
    def __init__(self, address, keys_file="", series="NumLock"):
        self.address = address
        try:
            with open(keys_file) as f:
                d = json.load(f)
            self.robot = d['robot']
            self.series = series
        except FileNotFoundError:
            keys_file = "keys.json"
            self.series = series
            self.register(keys_file)
        return

    def register(self, keys_file):
        r = requests.post(self.address+'get_sn/', data={
            'series': self.series,
        })
        self.robot = r.json()['SN']
        with open(keys_file, 'w+') as f:
            d = {
                "robot": self.robot,
                "series": self.series,
                "key": r.json()['key'],
                "jwt_key": r.json()['jwt_key']
            }
            json.dump(d, f)
        return r
