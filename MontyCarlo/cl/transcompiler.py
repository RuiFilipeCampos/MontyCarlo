


"""
asdas adasd
asd sad  # asdasd


<Sphere  xxx >



"""




class Token:
    def parse_symbol(self):
        pass

class Tag:
    def __init__(self):
        self.getting_name = True
        self.self_closing = None

        self.end = False
    

    def parse_symbol(self, symbol):
        if self.end:
            self.next = self.determine_next_token()
        
        if symbol == "/":
            self_closing = True
            return False

        if symbol == ">":
            self.end = True
            return False
        
        if self.getting_name:

            if symbol == " ":
                self.getting_name = False
                return False

            self.name += symbol

        if symbol == " ":
            return False

        
        self.getting_var_name = True

        if self.getting_var_name:
            

        if self.getting_value:
            pass

        
            
        


class PythonCode:
    def __init__(self):
        self.code = ""


    def __bool__(self):
        return True
    
    def parse_symbol(self, symbol):
        if symbol == "<":
            self.next = Tag()
            return True

        self.code += symbol
        return False


class SelfClosing:
    pass



def parse(code):

    tokens = []
    current_token = PythonCode()

    for letter in code:
        if current_token.parse_symbol(letter):
            current_token = current_token.next



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



            