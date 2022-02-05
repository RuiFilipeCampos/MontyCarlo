from django.db import models


class Directory(models.Model):
    name = models.CharField(max_length=500)
    project = models.ForeignKey(
        "projects.Project",
        on_delete=models.CASCADE,
        related_name='directories'
    )
    parent = models.ForeignKey(
        "projects.Directory",
        on_delete=models.CASCADE
    )


class CodeFile(models.Model):
    name = models.CharField(max_length=500)
    parent = models.ForeignKey(
        "projects.Directory",
        on_delete=models.CASCADE,
    )
    project = models.ForeignKey(
        "projects.Project",
        on_delete=models.CASCADE,
        related_name='code_files'
    )
    pythonX_code = models.TextField()
    python_code = models.TextField()

class Project(models.Model):
    name = models.CharField(max_length=500)
    # code_files = [ ]
    # directories = [ ]
