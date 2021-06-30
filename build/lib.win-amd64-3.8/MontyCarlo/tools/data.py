from numpy import array

def getAxis(text):
    """Get axis from list of lines of text."""
    
    xAxis, yAxis = [], []
    for line in text:
        XY = [float(number) for number in line.split()]
        xAxis += [XY[0]]
        yAxis += [XY[1]]
    xAxis, yAxis = map(array, (xAxis, yAxis))
    return xAxis, yAxis



def getTable(text):
    line0 = text[0]
    N = len(line0.split())
    table = []    
    for line in text:
         tableLine = array([float(number) for number in line.split()])
         table += [tableLine]
    table = array(table)
    return table


def get_bookmarked_text(path):
    """
    Reads EPICS file format and returns a dict with flags as keys
    and data as a list of strings.
    """
    with open(path, "r") as file:
        text = file.readlines()
        text = [line.strip('\n') for line in text]

        bookmarks = [0]
        
        for n, line in enumerate(text):
            if line == "                                                                       1":
                bookmarks += [n + 1]
                
        #gather all bookmarked text into a dict
        bookmark_ = bookmarks[0:-1]
        _ookmarks = bookmarks[1:]
    
        bookmarked_text = {}
        for i, j in zip(bookmark_, _ookmarks):
            line1, line2 = text[i], text[i+1]

            #on line 1
            Yi = float(line1[7:9])    #particle identifier
            Yo = float(line1[10:12])  #secondary particle designator
            Iflag = float(line1[31])  #interpolation flag
            
            #on line 2
            C  = float(line2[0:2])    #reaction descriptor
            I  = float(line2[2:5])    #reaction property
            S  = float(line2[5:8])    #reaction modifier
            X1 = float(line2[22:32])  #subshell designator
            
            flags = (Yi, C, S, X1, Yo, I)

            flags = tuple(map(int, flags))
            bookmarked_text[flags] = (Iflag, text[i+2:j-1])
            
    return bookmarked_text

