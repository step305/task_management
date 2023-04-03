# Generated by Django 3.0.5 on 2023-03-30 14:01

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('RobotMonitor', '0002_auto_20200603_1630'),
    ]

    operations = [
        migrations.AlterField(
            model_name='firmware',
            name='Archieve',
            field=models.FileField(upload_to='Firmware', verbose_name='Прошивка'),
        ),
        migrations.AlterField(
            model_name='firmware',
            name='Type',
            field=models.CharField(choices=[('VT', 'Прошивка видео для башни'), ('VM', 'Прошивка видео для базы'), ('ROS', 'Прошивка ROS')], max_length=3, verbose_name='Тип прошивки'),
        ),
    ]