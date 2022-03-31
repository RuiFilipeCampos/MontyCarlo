# pages/urls.py
from django.urls import path
from pages.views import IndexPageView

urlpatterns = [
    path("", IndexPageView.as_view(), name="index"),
]