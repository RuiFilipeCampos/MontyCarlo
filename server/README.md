


# Setting up for development.

Install [poetry]().

Run the following commands:

```
poetry shell
poetry install
python manage.py makemigrations
python manage.py migrate
```

# Running local dev server

```
python manage.py runserver
```

This defaults to port `8000` on localhost.
