

from collections import OrderedDict
from rest_framework import serializers
from materials.models import Material


class MaterialSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Material
        fields = "__all__"


    def to_representation(self, instance):
        result = super(MaterialSerializer, self).to_representation(instance)
        return OrderedDict([(key, result[key]) for key in result if result[key] ])
