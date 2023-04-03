from django.db import models
from django_countries.fields import CountryField

from django.utils.translation import gettext_lazy as _
from django.utils import timezone

class series(models.Model):
    Distribution_key = models.CharField('Наименование серии', max_length=255, unique=True)
    SN_amount = models.PositiveIntegerField('Количество роботов этой серии')
    SN_distributed = models.PositiveIntegerField('Количество выданных серийных номеров этой серии')
    update_time = models.DateTimeField('Дата изменения', auto_now = True)
    inv_status = models.BooleanField('Скрыть из общего списка', default = False)

    def __str__(self):
        return 'Серия: "{}", {}/{}'.format(self.Distribution_key, self.SN_distributed, self.SN_amount)


class robot(models.Model):
    serial_number = models.PositiveIntegerField('Серийный номер', unique = True)
    series = models.ForeignKey(series, on_delete = models.PROTECT, verbose_name = 'Серия', null = True)
    class statuses(models.TextChoices):
        Explotation = 'EX', _('Эксплуатация')
        OnRepair = 'OR', _('Ремонт')
        NeedRepair = 'NR', _('Требует ремонта')
        Broken = 'BR', _('Сломан')
        OnHold = 'OH', _('В ожидании')
    status = models.CharField(
        'Статус',
        max_length=2,
        choices=statuses.choices,
        default=statuses.OnHold,
    )
    ekey = models.CharField('Ключ', max_length=255)
    update_time = models.DateTimeField('Дата изменения', auto_now = True)
    inv_status = models.BooleanField('Скрыть из общего списка', default = False)

    def __str__(self):
        return 'Робот нормер {}, Статус {}'.format(self.serial_number, self.get_status_display())

class location(models.Model):
    name = models.CharField('Название', max_length=255)

    class types(models.TextChoices):
        Exhibition = 'EX', _('Выставка')
        Airport = 'AP', _('Аэропорт')
        Service = 'SE', _('Сервис')
        Storage = 'ST', _('Склад')
        Production = 'PR', _('Производство')
    status = models.CharField(
        'Статус',
        max_length=2,
        choices=types.choices
    )

    country = CountryField('Страна', blank=True, null=True)
    town = models.CharField('Город', max_length=255, blank=True, null=True)
    lon = models.FloatField('Долгота', blank=True, null=True)
    lat = models.FloatField('Широта', blank=True, null=True)
    date_in =  models.DateTimeField('Дата создания', auto_now_add = True)
    update_time = models.DateTimeField('Дата изменения', auto_now = True)
    inv_status = models.BooleanField('Скрыть из общего списка', default = False)
    def __str__(self):
        return self.name

class robot_belong(models.Model):
    robot = models.ForeignKey(robot, on_delete = models.PROTECT, verbose_name = 'Робот')
    company_name = models.CharField('Название компании', max_length=255)
    location = models.ForeignKey(location, on_delete = models.PROTECT, verbose_name = 'Место прибывания')

    class aims(models.TextChoices):
        Explotation = 'EX', _('Эксплуатация')
        Exhibition = 'EH', _('Показ')
        Repair = 'RP', _('Ремонт')
        Safekeeping = 'SK', _('Хранение')
        Utilization = 'UT', _('Утилизация')
    aim = models.CharField(
        'Цель прибывания',
        max_length=2,
        choices=aims.choices,
    )

    limitation = models.DateTimeField('Дата до которой отдан робот', blank=True, null=True)
    date_in = models.DateTimeField('Дата прибытия', blank=True, null=True)
    from_id = models.ForeignKey('self', on_delete = models.PROTECT, blank=True, null=True, verbose_name = 'Получено от')
    date_out =  models.DateTimeField('Дата отправки', blank=True, null=True)

    class statuses(models.TextChoices):
        Operational = 'OP', _('Рабочий')
        NeedRepair = 'NR', _('Требует ремонта')
        Broken = 'BR', _('Сломан')
    out_status = models.CharField(
        'Состояние в котором был отправлен',
        max_length=2,
        choices=statuses.choices,
        blank=True, null=True
    )

    update_time = models.DateTimeField('Дата изменения', auto_now = True)
    inv_status = models.BooleanField('Скрыть из общего списка', default = False)

    def __str__(self):
        return 'Компания {}, робот номер {}'.format(self.company_name, self.robot.serial_number)

class robot_session(models.Model):
    robot = models.ForeignKey(robot, on_delete = models.PROTECT, verbose_name = 'Робот')
    av_lon = models.FloatField('Долгота', blank=True, null=True)
    av_lat = models.FloatField('Широта', blank=True, null=True)
    volt_start = models.FloatField('Начальное напряжение', blank=True, null=True)
    volt_end = models.FloatField('Конечное напряжение', blank=True, null=True)
    GPS_status = models.IntegerField('Статус GPS', blank=True, null=True)
    track_file = models.FileField('Файл GPS трека', upload_to='tracks/', blank=True, null=True)
    media_file = models.FileField('Медиа файл', upload_to='media_from_robot/', blank=True, null=True)
    sensors_file = models.FileField('Файл данных сенсоров', upload_to='sensors/', blank=True, null=True)

    class ends(models.TextChoices):
        Continues = 'OL', _('На связи')
        Closed = 'CR', _('Корректное завершение работы')
        Connection_lost = 'DC', _('Потеря связи')
    out_status = models.CharField(
        'Тип завершения сессии',
        max_length=2,
        choices=ends.choices
    )

    date_time_start = models.DateTimeField('Начало записи', auto_now_add = True)
    date_time_end = models.DateTimeField('Конец записи', auto_now = True)
    update_time = models.DateTimeField('Дата изменения', auto_now = True)
    inv_status = models.BooleanField('Скрыть из общего списка', default = False)

    def __str__(self):
        return 'Робот номер {} Записано: '.format(self.robot.serial_number)+self.date_time_start.strftime("%Y-%m-%d %H:%M:%S")

class sensor(models.Model):
    robot = models.ForeignKey(robot, on_delete = models.PROTECT, verbose_name = 'Робот')
    session = models.ForeignKey(robot_session, on_delete = models.PROTECT, verbose_name = 'Сессия')
    sensor = models.CharField('Название сенсора', max_length=50)
    node = models.BooleanField('Программа', default = False)
    status = models.BooleanField('Рабочий', default = False)
    correctness = models.BooleanField('Корректные данные', default = False)

    class types(models.TextChoices):
        Auto = 'AT', _('Автоматическая запись')
        Demand = 'DM', _('Запись по требованию')
    out_status = models.CharField(
        'Тип записи',
        max_length=2,
        choices=types.choices
    )

    date_time_in = models.DateTimeField('Время прихода информации', auto_now_add = True)
    update_time = models.DateTimeField('Дата изменения', auto_now = True)
    inv_status = models.BooleanField('Скрыть из общего списка', default = False)

    def __str__(self):
        st = 'рабочий' if self.status else 'не рабочий'
        return 'Робот нормер {}, Сенсор {}, Статус {}'.format(self.robot.serial_number, self.sensor, st)


class firmware(models.Model):
    Version_number = models.CharField('Номер версии', max_length=50)
    class types(models.TextChoices):
        Video_top = 'VT', _('Прошивка видео для башни')
        Video_main = 'VM', _('Прошивка видео для базы')
        ROS = 'ROS', _('Прошивка ROS')
    Type = models.CharField(
        'Тип прошивки',
        max_length=3,
        choices=types.choices
    )
    Archieve = models.FileField('Прошивка', upload_to='Firmware')
    update_time = models.DateTimeField('Дата изменения', auto_now = True)
    inv_status = models.BooleanField('Скрыть из общего списка', default = False)
    def __str__(self):
        return '{} версии {}'.format(self.get_Type_display(), self.Version_number)
