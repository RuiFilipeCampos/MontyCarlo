


"""
asdas adasd
asd sad  # asdasd


<Sphere>



"""


def to_python(code):

    tokens = []
    inside_tag = False

    writing_name = False
    writing_func_name = False

    for letter in code:
        if inside_tag:


            if writing_func_name:
                if letter == " ":
                    writing_func_name = False
                    continue

                if letter == ">":
                    writing_func_name = False
                    tokens.append(">")
                    continue

                tokens[-1] += letter
                continue



            if letter == " ":
                continue
            
            if letter == "=":
                tokens.append("=")
                tokens.append("")
                continue

            if letter == "/":
                tokens.append("/")
                continue

            if letter == ">":
                tokens.append(">")
                inside_tag = False
                continue



            



        if letter == "<":
            tokens.append("<")
            tokens.append("")
            inside_tag = True
            writing_func_name = True


    return tokens


if __name__ == "__main__":

    print(to_python(r"""
    def f(x):
        return x**2
    
    <Hello> <Hello>
    
    
    
    """))



            