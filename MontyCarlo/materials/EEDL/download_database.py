import requests

def construct_url(n: int) -> str:
    if n < 10:
        return f"https://www-nds.iaea.org/epics/ENDL2017/EEDL.ELEMENTS/ZA00{n}000"

    if n == 100:
        return "https://www-nds.iaea.org/epics/ENDL2017/EEDL.ELEMENTS/ZA100000"

    return  f"https://www-nds.iaea.org/epics/ENDL2017/EEDL.ELEMENTS/ZA0{n}000"

for N in range(1, 101):
    print(f"Downloading EEDL {N}")

    url =  construct_url(N)
    response = requests.get(url)
    content = response.content

    with open(f"{N}.txt", "wb") as file:
        file.write(content)

