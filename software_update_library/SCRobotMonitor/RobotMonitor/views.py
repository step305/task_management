from django.shortcuts import render, get_object_or_404, Http404
from django.views.decorators.csrf import csrf_exempt
from django.http import JsonResponse, FileResponse
from django.db.models import Max

import jwt, base64, json, random, string

from .models import sensor, robot, robot_session, series, firmware
from .pymodules.SC_encryption import SC_encryption


with open('RobotMonitor/pymodules/keys.json') as f:
    d = json.load(f)
encryption = SC_encryption(d['jwt_key'])
encryption.set_key(encryption.decode_key(d["key"]))

rand_str = lambda n: ''.join([random.choice(string.ascii_lowercase) for i in range(n)])

# Create your views here.
def getSN(request):
    try:
        series_name = request.POST['series']
        series_object = series.objects.get(Distribution_key=series_name)
        if series_object.SN_amount>series_object.SN_distributed:
            try:
                new_number = robot.objects.aggregate(Max('serial_number'))['serial_number__max']+1

                series_object.SN_distributed+=1
            except:
                new_number = 7
                new_key = encryption.gen_new_key()
                series_object.SN_distributed+=1
        else:
            return JsonResponse({'Result':'No free SNs'})
    except:
        return JsonResponse({'Result':'Series not found'})
    else:
        new_key = encryption.gen_new_key()
        new_robot = robot(
            serial_number = new_number,
            series = series_object,
            status = 'OH',
            ekey = encryption.encode_key(new_key)
        )
        new_robot.save()
        series_object.save()
        return JsonResponse({'Result':'Added', 'SN': new_number, 'key':encryption.encode_key(new_key), 'jwt_key':encryption.get_jwt_key()})

def startSession(request):
    try:
        token = request.POST['data']
        decoded = encryption.decode_jwt(token)
        robot_id = decoded['data']
        try:
            recieved_from = robot.objects.get(serial_number=robot_id)
            encryption.set_key(encryption.decode_key(recieved_from.ekey))
            try:
                older_session = robot_session.objects.filter(robot=recieved_from).filter(out_status='OL').order_by('-update_time')
                if len(older_session)>0:
                    older_session = older_session[0]
            except (KeyError):
                older_session = False
        except (KeyError):
            return JsonResponse({'Result':'Robot not found'})
    except (KeyError):
        # Redisplay the question voting form.
        return JsonResponse({'Result':'Error wrong request'})
    else:
        if older_session:
            older_session.out_status='DC'
            older_session.save()
        s = robot_session(
            robot = recieved_from,
            out_status = 'OL'
        )
        s.save()
        new_key = encryption.gen_new_key()
        recieved_from.ekey = encryption.encode_key(new_key)
        recieved_from.save()
        data = encryption.code_data({'key':encryption.encode_key(new_key)})
        encryption.set_key(new_key)
        return JsonResponse({'Result':'Session started', 'data': data})

def addSensor(request):
    try:
        eej = request.POST['data']
        robot_id = request.POST['robot']
        try:
            coded = encryption.decode_jwt(eej)
            try:
                recieved_from = robot.objects.get(serial_number=robot_id)
                encryption.set_key(recieved_from.ekey)
            except (KeyError):
                return JsonResponse({'Result':'Robot not found'})
            data = encryption.decode_data(coded['data'])
            sensor_name = data['sensor']
            node = data['node']
            valid = data['topic']
            valid_data = data['data']
        except:
            print('Wrong sig')
        try:
            current_session = robot_session.objects.filter(robot=recieved_from).order_by('-update_time')[0]
        except (KeyError):
            return JsonResponse({'Result':'No sessions for this robot'})
    except (KeyError):
        return JsonResponse({'Result':'Error wrong request'})
    else:
        s = sensor(
            robot = recieved_from,
            session = current_session,
            sensor = sensor_name,
            node = node,
            status = valid,
            correctness = valid_data,
            out_status = 'AT',
        )
        s.save()
        return JsonResponse({'Result':'Sensor saved'})

def addSession(request):
    try:
        eej = request.POST['data']
        robot_id = request.POST['robot']
        try:
            coded = encryption.decode_jwt(eej)
            try:
                recieved_from = robot.objects.get(serial_number=robot_id)
                encryption.set_key(recieved_from.ekey)
            except (KeyError):
                return JsonResponse({'Result':'Robot not found'})
            data = encryption.decode_data(coded['data'])
            av_lon = data['av_lon']
            av_lat = data['av_lat']
            volt_start = data['volt_start']
            volt_end = data['volt_end']
            GPS_status = data['GPS_status']
        except:
            print('Wrong sig')
        try:
            current_session = robot_session.objects.filter(robot=recieved_from).order_by('-update_time')[0]
        except (KeyError):
            return JsonResponse({'Result':'No sessions for this robot'})
    except (KeyError):
        return JsonResponse({'Result':'Error wrong request'})
    else:
        if av_lon != 0:
            current_session.av_lon = av_lon
        if av_lat != 0:
            current_session.av_lat = av_lat
        if volt_start != -1:
            current_session.volt_start = volt_start
        if volt_end != -1:
            current_session.volt_end = volt_end
        if GPS_status != -1:
            current_session.GPS_status = GPS_status
        current_session.save()
        return JsonResponse({'Result':'Session data saved'})

def closeSession(request):
    try:
        token = request.POST['data']
        decoded = encryption.decode_jwt(token)
        robot_id = decoded['data']
        try:
            recieved_from = robot.objects.get(serial_number=robot_id)
            encryption.set_key(recieved_from.ekey)
        except (KeyError):
            return JsonResponse({'Result':'Robot not found'})
        try:
            current_session = robot_session.objects.filter(robot=recieved_from).order_by('-update_time')[0]
        except (KeyError):
            return JsonResponse({'Result':'No sessions for this robot'})
    except (KeyError):
        return JsonResponse({'Result':'Error wrong request'})
    else:
        if current_session.out_status == 'OL':
            current_session.out_status = 'CR'
            current_session.save()
            return JsonResponse({'Result':'Session closed'})
        else:
            return JsonResponse({'Result':'Session already closed'})

def uploadSessions(request):
    try:
        token = request.POST['data']
        file = request.FILES['media']
        decoded = encryption.decode_jwt(token)
        robot_id = decoded['data']
        try:
            recieved_from = robot.objects.get(serial_number=robot_id)
            encryption.set_key(recieved_from.ekey)
        except (KeyError):
            return JsonResponse({'Result':'Robot not found'})
        try:
            current_session = robot_session.objects.filter(robot=recieved_from).order_by('-update_time')[0]
        except (KeyError):
            return JsonResponse({'Result':'No sessions for this robot'})
    except (KeyError):
        return JsonResponse({'Result':'Error wrong request'})
    else:
        if current_session.out_status == 'OL':
            if 'bag' in request.FILES:
                current_session.sensors_file = request.FILES['bag']
            if 'media' in request.FILES:
                current_session.media_file = request.FILES['media']
            if 'GPS' in request.FILES:
                current_session.track_file = request.FILES['GPS']
            current_session.save()
            return JsonResponse({'Result':'Files saved'})
        else:
            return JsonResponse({'Result':'Session already closed'})

def getFirmwareVersion(request):
    try:
        token = request.POST['data']
        decoded = encryption.decode_jwt(token)
        robot_id = decoded['data']
        try:
            recieved_from = robot.objects.get(serial_number=robot_id)
            encryption.set_key(recieved_from.ekey)
        except (KeyError):
            return JsonResponse({'Result':'Robot not found'})
        try:
            firmwareType = decoded['type']
        except (KeyError):
            return JsonResponse({'Result':'Error wrong request1'})
    except (KeyError):
        return JsonResponse({'Result':'Error wrong request2'})
    else:
        try:
            if firmwareType == 'ROS':
                firmwareObj = firmware.objects.filter(Type='ROS').order_by('-update_time')[0]
            elif firmwareType == 'TowerVideo':
                firmwareObj = firmware.objects.filter(Type='VT').order_by('-update_time')[0]
            elif firmwareType == 'MainVideo':
                firmwareObj = firmware.objects.filter(Type='VM').order_by('-update_time')[0]
            else:
                return JsonResponse({'Result':'Error wrong request3'})
        except (KeyError):
            return JsonResponse({'Result':'Error wrong request4'})
        else:
            return JsonResponse({'Result':firmwareObj.Version_number})

def getFirmware(request):
    try:
        token = request.POST['data']
        decoded = encryption.decode_jwt(token)
        robot_id = decoded['data']
        try:
            recieved_from = robot.objects.get(serial_number=robot_id)
            encryption.set_key(recieved_from.ekey)
        except (KeyError):
            return JsonResponse({'Result':'Robot not found'})
        try:
            firmwareType = decoded['type']
            versionNumber = decoded['version']
        except (KeyError):
            return JsonResponse({'Result':'Error wrong request'})
    except:
        return JsonResponse({'Result':'Error wrong request'})
    else:
        try:
            if firmwareType == 'ROS':
                firmwareObj = firmware.objects.filter(Type='ROS').filter(Version_number = versionNumber)[0]
            elif firmwareType == 'TowerVideo':
                firmwareObj = firmware.objects.filter(Type='VT').filter(Version_number = versionNumber)[0]
            elif firmwareType == 'MainVideo':
                firmwareObj = firmware.objects.filter(Type='VM').filter(Version_number = versionNumber)[0]
            else:
                return JsonResponse({'Result':'Error wrong request'})
        except (KeyError):
            return JsonResponse({'Result':'Error wrong request'})
        else:
            return FileResponse(open(firmwareObj.Archieve.name, 'rb'))
