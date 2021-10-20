


"""
asdas adasd
asd sad  # asdasd


<Sphere  xxx >



"""



class Variable:
    def __init__(self, symbol):
        self.name = symbol
        self.value = ""



DETERMINING_CLOSURE = 0
GETTING_NAME = 1
SEARCHING = 2
GETTING_VAR_NAME = 3
GETTING_VAR_VALUE = 4

ENDING_TOKEN = -1
class Tag:


    def __init__(self):

        self.state = DETERMINING_CLOSURE
        self.vars = []

        self.name = ""


        self.type = ""

    def __repr__(self):
        if self.type == "self-closing":
            return f"<{self.name} />"
        if self.type == "start":
            return f"<{self.name}>"
        if self.type == "end":
            return f"</{self.name}>"
        
        return f"<{self.name}>"

    def to_code(self):
        if self.vars:
            arguments = ""
            for var in self.vars:
                arguments += f"{var.name}={var.value},"

            return f"{self.name}({arguments})"
        return f"{self.name}()"

    def determine_next_token(self, symbol):
        if symbol == "<":
            return Tag()
    
        return PythonCode(symbol)

    def parse_symbol(self, symbol):

        if self.state == ENDING_TOKEN:
            self.next = self.determine_next_token(symbol)
            return True

        if self.state == DETERMINING_CLOSURE:
            if symbol == "/":
                self.type = "end"
                self.getting_name = True
            else:
                self.type = "start"
                self.name += symbol

            self.state = GETTING_NAME
            return False
        

        # CHECK FOR END
        if symbol == "/":
            self.type = "self-closing"
            return False

        if symbol == ">":
            self.state = ENDING_TOKEN
            return False
        
        if self.state == GETTING_NAME:

            if symbol == " ":
                self.state = SEARCHING
                return False

            self.name += symbol
            return False

        if self.state == SEARCHING:
            if symbol == " ":
                return False
            
            self.state = GETTING_VAR_NAME
            self.current_variable = Variable(symbol)
            return False

        

        if self.state == GETTING_VAR_NAME:
            if symbol == " ":
                raise RuntimeError("")

            if symbol == "=":
                self.state = GETTING_VAR_VALUE
                return False
            
            self.current_variable.name += symbol
            return False
            

        if self.state == GETTING_VAR_VALUE:
            if symbol == " ":
                self.state = SEARCHING
                self.vars.append(self.current_variable)
                return False 

            self.current_variable.value += symbol
            return False

        
            
        


class PythonCode:
    def __init__(self, symbol):
        self.code = symbol

    def to_code(self):
        return self.code
    
    def __bool__(self):
        return True

    def __repr__(self):
        return "<PythonCode>"
    
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
    current_token = PythonCode("")
    tokens.append(current_token)

    for letter in code:
        print(current_token, letter)
        if current_token.parse_symbol(letter):
            

            current_token = current_token.next
            tokens.append(current_token)
    return tokens


# are these names even correct? xdd
# gotta check

def to_code(tag):
    code = """
    def {tag}
    """




if __name__ == "__main__":
    code = r"""
    def f(x):
        return x**2

    def g(x):
        return (
            <Hummm />
)
    
    <Hello/>

    <World><HI/>
    </World>

    
    
    
    """

    print("Input:")
    print(code)
    print("")

    print("Output:")
    TOKENS = parse(code)
    print(TOKENS)

    print("")

    for token in TOKENS:
        print(token.to_code())










