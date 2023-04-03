from django.contrib import admin

from .models import robot, location, robot_belong, robot_session, sensor, series, firmware

admin.site.register(series)
admin.site.register(robot)
admin.site.register(location)
admin.site.register(robot_belong)
admin.site.register(robot_session)
admin.site.register(sensor)
admin.site.register(firmware)
