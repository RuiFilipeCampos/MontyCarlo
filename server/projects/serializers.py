



from collections import OrderedDict
from rest_framework import serializers
from projects.models import Project
from projects.models import CodeFile



class ProjectSerializer(serializers.HyperlinkedModelSerializer):
    code_files  = serializers.RelatedField(many=True, read_only=True)
    directories = serializers.RelatedField(many=True, read_only=True)
    class Meta:
        model = Project
        fields = ['id', 'name', 'code_files', 'directories']

