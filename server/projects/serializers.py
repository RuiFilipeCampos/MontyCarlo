



from collections import OrderedDict
from rest_framework import serializers
from projects.models import Project
from projects.models import CodeFile
from projects.models import Directory


from rest_framework_recursive.fields import RecursiveField



class CodeFileSerializer(serializers.ModelSerializer):
    class Meta:
        model = CodeFile
        fields = '__all__'

class DirectorySerializer(serializers.ModelSerializer):
    # code_files  = CodeFileSerializer(many=True, read_only=True)
    directories = serializers.PrimaryKeyRelatedField(
        read_only=True,
        many=True
    )
    code_files  = CodeFileSerializer(many=True, read_only=True)

    class Meta:
        model = Directory
        fields = '__all__'

class ProjectSerializer(serializers.ModelSerializer):
    code_files  = CodeFileSerializer(many=True, read_only=True)
    directories = DirectorySerializer(many=True, read_only=True)
    
    
    
    # serializers.PrimaryKeyRelatedField(many=True, read_only=True)
    class Meta:
        model = Project
        fields = '__all__'

    
    def to_representation(self, instance):
        representation = super(ProjectSerializer, self).to_representation(
            instance
        )


        new_directories = []
        for folder in representation['directories']:
            if folder['parent'] == None:
                new_directories.append(folder['id'])
        representation['directories'] = new_directories

        new_files = []
        for file in representation['code_files']:
            if file['parent'] == None:
                new_directories.append(file)
        representation['code_files'] = new_files
        return representation