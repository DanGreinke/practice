from fastapi import FastAPI
import redis

app = FastAPI()

r = redis.Redis(host='redis', port=6379, db=0)

import debugpy
debugpy.listen(("0.0.0.0", 5678))

@app.get("/")
def read_root():
    return {"Hello": "World123"}

@app.get("/hits")
def read_root():
    r.incr("hits")
    return {"Number of hits": r.get("hits")}