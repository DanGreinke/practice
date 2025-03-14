# Docker Tutorial
Followed from youtube video, with Patrick Loeber
https://www.youtube.com/watch?v=6OxqiEeCvMI 
## Build and Run Image

To build the image, enter:

> $ docker build -t fastapi-image .

To run image (without docker compose):

> $ docker run --name fastapi-container -p 80:80 -d -v $(pwd):/code fastapi-image

Details on run cmd
* -p 80:80 maps port 80 in the container to the host machine's port 80
* -d for detached mode so we can use the terminal after the container starts running
* -v $(pwd):/code maps current working dir in the host machine to /code dir in container, so we can modify code while container is running
* fastapi-image specifies the image used to build the container

## Use VSCode in Container
* Install Docker and Dev Container extensions
* Can attach to open container from bottom left corner button ><
* Intall python extension on VSCode instance inside the container
* Changes made inside container are propagated outside to the host machine

## Stop & Remove Container

> $ docker stop fastapi-container

> $ docker rm fastapi-container

## Using docker compose to specify params and simplify container start

* See docker-compose.yml for params

To start container, enter:
> $ docker-compose up

To remove container, enter:
> $ docker-compose down

## Added Redis to docker-compose
* Service added
* Specified dependency on redis in app
* Redis will be our database
* Added decorator for hits, to increment when we refresh the page.

To Build
> $ docker-compose up --build -d

## Add debugpy
* Add port 5678 to docker compose
* Add debugpy to reqs
* Add python code to main for debugpy
* Run debugger, create json config, Python Debugger > Remote Attach
* Can now run python in debug mode and pause on breakpoints


End