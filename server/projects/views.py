

from rest_framework import viewsets

from projects.models import Project, Directory, CodeFile
from projects.serializers import ProjectSerializer
from projects.serializers import DirectorySerializer
from projects.serializers import CodeFileSerializer


class ProjectViewSets(viewsets.ModelViewSet):
    queryset = ProjectSerializer.Meta.model.objects.all()
    serializer_class = ProjectSerializer


class DirectoryViewSets(viewsets.ModelViewSet):
    queryset = DirectorySerializer.Meta.model.objects.all()
    serializer_class = DirectorySerializer

    def get_queryset(self):
        project = Project.objects.get(pk=int(self.kwargs['project_pk']))
        return DirectorySerializer.Meta.model.objects.filter(project=project)

class CodeFileViewSets(viewsets.ModelViewSet):
    queryset = CodeFileSerializer.Meta.model.objects.all()
    serializer_class = CodeFileSerializer

    def get_queryset(self):
        project = Project.objects.get(pk=int(self.kwargs['project_pk']))
        return CodeFileSerializer.Meta.model.objects.filter(project = project)