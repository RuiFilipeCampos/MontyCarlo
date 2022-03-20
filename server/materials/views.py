from django.shortcuts import render

from rest_framework import viewsets


from materials.serializers import MaterialSerializer


class MaterialViewSets(viewsets.ModelViewSet):
    queryset = MaterialSerializer.Meta.model.objects.all()
    serializer_class = MaterialSerializer