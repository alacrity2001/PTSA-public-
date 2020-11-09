#######################################################################
#######################################################################
#####				MATRIX PROFILE PLOT FUNCTIONS				 ######
#######################################################################
#######################################################################


import numpy as np
import matplotlib.pyplot as plt
import plotly
import plotly.graph_objs as go


def plot(ts, title="Matrix Profile",type_ts='mp',figsize=(20,10),fontlevel=2, point_size=100):
    """
        Return a tuple(fig,ax) where figure is the figure to show or store
		and ax is the axes of the matplotlib plot.
    """
    if type_ts not in ['mp','ip']:
    	raise ValueError("type_ts must be in ['mp','ip']")

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.set_title(title,fontsize=fontlevel*8)
    ax.plot(ts,color="#166aca",alpha=0.7)
    if type_ts == 'mp':
        ax.scatter(np.argmax(ts),max(ts),color="#1bd048",s=point_size,label="top discord")
        ax.scatter(np.argmin(ts),min(ts),color="#c91818",s=point_size,label="top motif")
        ax.legend(fontsize=fontlevel*8)
        ax.set_ylabel("Minimum Distance",fontsize=fontlevel*6)
    elif type_ts == 'ip':
        ax.set_ylabel("Index position nearest neighbor",fontsize=fontlevel*6)

    ax.set_xlim(0,len(ts))
    ax.set_xlabel("Index",fontsize=fontlevel*6)
    ax.set_ylim(0,np.max(ts)+2)
    return fig,ax

def i_plot(ts, title="Matrix Profile",type_ts='mp',figsize=(950,500),fontlevel=2, point_size=10):
    """
        Return a plotly figure that allows to create interactive plot
    """
    if type_ts not in ['mp','ip']:
        raise ValueError("type_ts must be in ['mp','ip']")


    data = []
    ########## Trace Matrix Profile ############
    data.append(
        go.Scattergl(
            x = np.arange(0,len(ts)),
            y = np.array(ts),
            opacity = 0.7,
            line = dict(color = "#166aca"),
            name=title
        )
    )
    if type_ts == 'mp':
        ylabel = 'Minimum Distance'
        ######### Trace Top Discord #############
        data.append(
            go.Scattergl(
                x = [np.argmax(ts)],
                y = [max(ts)],
                mode = 'markers',
                marker = dict(color = "#1bd048", size = point_size),
                name='Top Discord'
            )
        )
        ######### Trace Top Motif #############
        data.append(
            go.Scattergl(
                x = [np.argmin(ts)],
                y = [min(ts)],
                mode = 'markers',
                marker = dict(color = "#c91818", size = point_size),
                name='Top Motif'
            )
        )
    elif type_ts == 'ip':
        ylabel = "Index position nearest neighbor"
    ######### Layout #############
    layout = dict(
        legend=dict(orientation="h"),
        height=figsize[1],
        width=figsize[0],
        title=title,
        titlefont=dict(size=fontlevel*9),
        xaxis=dict(
            title='Index',
            titlefont=dict(
                size=fontlevel*6,
            )
        ),
        yaxis=dict(
            title=ylabel,
            titlefont=dict(
                size=fontlevel*6,
            )
        ),
    )

    return dict(data=data,layout=layout)
