import jwt

encodedJwt = jwt.encode({"some": "payload"}, "secret", algorithm="HS256")
print(encodedJwt, end="")
