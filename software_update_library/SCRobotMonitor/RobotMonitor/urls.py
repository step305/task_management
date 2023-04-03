from django.urls import path

from . import views

app_name = 'RobotMonitor'
urlpatterns = [
    # path('', views.IndexView.as_view(), name='index'),
    path('report_sensor/', views.addSensor, name='Загрузка данных сенсора'),
    path('report_session/', views.addSession, name='Загрузка данных сессии'),
    path('start_session/', views.startSession, name='Создание новой сессии'),
    path('close_session/', views.closeSession, name='Закрытие сессии'),
    path('session_upload/', views.uploadSessions, name='Загрузка файлов сессии'),
    path('get_sn/', views.getSN, name='Получение серийного номера'),
    path('get_latest_firmware_version/', views.getFirmwareVersion, name='Получение номера последней версии ПО'),
    path('get_latest_firmware/', views.getFirmware, name='Получение ПО'),
]
