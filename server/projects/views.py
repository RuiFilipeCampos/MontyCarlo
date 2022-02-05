

from rest_framework import viewsets


from projects.serializers import ProjectSerializer


class ProjectViewSets(viewsets.ModelViewSet):
    queryset = ProjectSerializer.Meta.model.objects.all()
    serializer_class = ProjectSerializer