from .htmlcreator import HTMLDocument
import numpy as np
from pandas import *




class MaterialLogger(HTMLDocument):
    def __init__(self, formula, density, name = "Untitled", 
                 C1 = 0.1, C2 = 0.1,
                 Wcr = 10e3, Wcc = 100e3):
        
        super(MaterialLogger, self).__init__()

        self.set_title(f"{name} - MyCo Material Output")

        self.add_header(f"Logging Material - {name}", level = "h1")

        self.add_header("Input Information")
        self.add_paragraph(f"Formula = {formula}")
        self.add_paragraph(f"Density = {density} g/cm^3")
        self.add_paragraph(f"C1 = {C1}")
        self.add_paragraph(f"C2 = {C2}")
        self.add_paragraph(f"wcc = {Wcc}eV")
        self.add_paragraph(f"Wcr = {Wcr}eV")
        

    def new_plot(self):

        import plotly.graph_objects as go
        return go.Figure()       
        

    def add_to_plot(self, fig, x, y, label = "no_label"):
        import plotly.graph_objects as go
        import numpy as np
        x = np.array(x)
        y = np.array(y)

        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            name = label
        ))



    def finish_plot(self, fig,  title = "Untitled", xlabel = "x-axis", ylabel = "y-axis", logscale = True, ylogscale = False, xlogscale = False):

        fig.update_layout(
            title=title,
            xaxis_title=xlabel,
            yaxis_title=ylabel,
            font=dict(
                family="Courier New, monospace",
                size=18,
                color="RebeccaPurple"
            )
        )


        if logscale:
            fig.update_xaxes(type="log")
            fig.update_yaxes(type="log")
        if ylogscale:
            fig.update_yaxes(type="log")
        if xlogscale:
            fig.update_xaxes(type="log")


        self.add_plotly_figure(fig)


    def add_plot(self, x, y, title = "Untitled", xlabel = "x-axis", ylabel = "y-axis", logscale = True, xlogscale = False, ylogscale = False):
        import numpy as np
        x = np.array(x)
        y = np.array(y)

        import plotly.graph_objects as go
        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=x,
            y=y,
        ))


        fig.update_layout(
            title=title,
            xaxis_title=xlabel,
            yaxis_title=ylabel,
            font=dict(
                family="Courier New, monospace",
                size=18,
                color="RebeccaPurple"
            )
        )


        if logscale:
            fig.update_xaxes(type="log")
            fig.update_yaxes(type="log")
        if xlogscale:
            fig.update_xaxes(type="log")
        if ylogscale:
            fig.update_yaxes(type="log")

        self.add_plotly_figure(fig)


    def add_attribute(self, name, value):
        self.add_paragraph(f"{name} = {value}")
    
    def log_table(self, xAxis, yAxis):
        df = DataFrame(data = (xAxis, yAxis))
        df.style.hide_index()
        self.add_table(df)
        

