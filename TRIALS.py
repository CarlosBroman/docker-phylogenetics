import plotly.graph_objects as go

values = [(13, 562, 58, 61), (6, 120, 0, 852), (459, 118, 16, 128), (68, 798, 68, 45), (92, 957, 81, 19), (50, 99, 1105, 3)]
positions = [1, 2, 3, 4, 5, 6]

# Create traces for each category
traces = []
categories = ['A', 'T', 'G', 'C']

for i, category in enumerate(categories):
    trace = go.Scatter(
        x=positions,
        y=[value[i] for value in values],
        mode='lines+markers',
        name=category
    )
    traces.append(trace)

# Create the layout
layout = go.Layout(
    title='Values for Different Categories',
    xaxis=dict(title='Position'),
    yaxis=dict(title='Value')
)

# Create the figure
fig = go.Figure(data=traces, layout=layout)

# Show the figure
fig.show()




