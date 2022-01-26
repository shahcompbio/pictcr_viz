FROM python:3.9.1
ADD . /python-flask
WORKDIR /python-flask
RUN pip install --default-timeout=100 -r requirements.txt
RUN export FLASK_APP=server.py
ENV FLASK_APP=/python-flask/server.py

CMD ["python3", "-m", "flask", "run", "--host=0.0.0.0"]
#ENTRYPOINT FLASK_APP=/python-flask/server.py flask run
#ENTRYPOINT [ "python" ]
#CMD [ "-u", "server.py" ]
