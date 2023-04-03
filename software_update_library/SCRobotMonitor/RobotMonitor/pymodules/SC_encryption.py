import jwt
import base64
import json
from cryptography.fernet import Fernet
# from Crypto.Cipher import DES


class SC_encryption():
    def __init__(self, jwt_key='', key=''):
        self.jwt_key = jwt_key
        self.key = key

    def pad(self, text):
        text = bytes(text, 'utf-8')
        while len(text) % 8 != 0:
            text += b' '
        return text

    def encode_jwt(self, data):
        jwt_encoded = jwt.encode(data, self.jwt_key, algorithm='HS512')
        return jwt_encoded

    def decode_jwt(self, token):
        coded = jwt.decode(token, self.jwt_key, algorithms='HS512')
        return coded

    def code_data(self, ob):
        des = Fernet(self.key)
        line = json.dumps(ob)
        spaced = self.pad(line)
        coded = des.encrypt(spaced)
        b64 = base64.encodebytes(coded).decode('ascii')
        return b64

    def decode_data(self, encoded):
        b = base64.b64decode(encoded)
        des = Fernet(self.key)
        line = des.decrypt(b)
        ob = json.loads(line)
        return ob

    def set_jwt_key(self, jwt_key):
        self.jwt_key = jwt_key

    def get_jwt_key(self):
        return self.jwt_key

    def set_key(self, key):
        self.key = key

    def gen_new_key(self):
        return Fernet.generate_key()

    def encode_key(self, key):
        return base64.urlsafe_b64encode(key).decode()

    def decode_key(self, key):
        return base64.urlsafe_b64decode(key.encode())
