import base64
import io
import dash
import pandas as pd
from Bio import SeqIO
from dash import Input, Output, State, dcc, html
import dash_bootstrap_components as dbc
import plotly.graph_objs as go

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SPACELAB])
server = app.server
app.title = 'DNA Sequence Analyser'

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H1("DNA Sequence Analyser")
        ], width=8)
    ]),
    dbc.Row([
        dbc.Col(dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Files')
            ]),
            style={
                'width': '60%',
                'height': '80px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            multiple=False)),
        dbc.Col([
            dcc.Store(id='data-store'),
        ], width=5)
    ]),
    dbc.Row([
        dbc.Col([
            html.Div(id='number-of-sequences', children=[], style={'font-size': '20px', 'font-weight': 'bold'}),
        ]),
    ]),
    dbc.Row([
        dbc.Col([
            dcc.Dropdown(
                id='dropdown',
                options=[],
                value='',
                placeholder='Select a sequence to analyse'),
        ], width=20)
    ]),
    dbc.Row([
        dbc.Col([
            dcc.Graph(
                id='seq-barplot'),
        ], width=6)
    ]),
    dbc.Row([
        dbc.Col([
            html.Div(id='output-sequence', children=[], style={'font-size': '20px', 'font-weight': 'bold'}),
        ]),
        dbc.Col([
            html.Div([
                dcc.Input(id='search-input', type='text', placeholder='Enter query subsequence'),
                html.Button('Search', id='search-button')
            ]),
        ])
    ]),
    dbc.Row([
        dbc.Col([
            html.Div(id='high-GC-subsequence', children=[], style={'font-size': '20px', 'font-weight': 'bold'}),
        ]),
        dbc.Col([
            html.Div(id='search-output', children=[], style={'font-size': '20px', 'font-weight': 'bold'})
        ])
    ])
])


@app.callback(
    Output('data-store', 'data'),
    Input('upload-data', 'contents'),
    State('upload-data', 'filename')
)
def store_data(contents, filename):
    if contents is not None:
        # Converts the uploaded CSV file to Pandas DataFrame
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        try:
            if 'csv' in filename:
                # Assume that the user uploaded a CSV file
                df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))

            elif 'fasta' or 'fa' in filename:
                records = SeqIO.parse(filename, "fasta")
                data = []
                for record in records:
                    data.append({"ID": record.id, "Sequence": str(record.seq)})
                df = pd.DataFrame(data)

            else:
                # Assume that the user uploaded an Excel file
                df = pd.read_excel(io.BytesIO(decoded))
        except Exception as e:
            print(e)
            return None

        # Store the DataFrame in dcc.Store component
        return df.to_dict('records')
    else:
        return None


@app.callback(
    Output('number-of-sequences', 'children'),
    Input('data-store', 'data'),
)
def no_seq_in_file(data):
    if data is not None:
        df = pd.DataFrame.from_dict(data)
        num_sequences = df.iloc[:, 0].nunique()
        output = "Number of sequences in uploaded file: {}".format(num_sequences)
        return output
    else:
        return None


@app.callback(
    Output('dropdown', 'options'),
    Input('data-store', 'data')
)
def update_dropdown_options(data):
    if data is not None:
        # Convert the stored data from dictionary to Pandas DataFrame
        df = pd.DataFrame.from_dict(data)
        # Generate options for the dropdown based on the sequence number and corresponding value from the second column
        options = [{'label': f'{seq} - {val}', 'value': val} for seq, val in zip(df.iloc[:, 0], df.iloc[:, 1])]
        # Set the default value of the dropdown to the first value in the options list
        default_value = options[0]['value']
        return options
    else:
        return []


@app.callback(
    Output('output-sequence', 'children'),
    Input('dropdown', 'value'),
    State('data-store', 'data')
)
def display_sequence_length(value, data):
    if data is not None and value is not None:
        # Convert the stored data from dictionary to Pandas DataFrame
        df = pd.DataFrame.from_dict(data)
        # Filter the DataFrame to get the row with the selected value
        mask = df.iloc[:, 1] == value
        filtered_df = df[mask]
        # Get the corresponding sequence from the second column
        sequence = filtered_df.iloc[0, 1]
        # Return the sequence with a message
        message = "Length of selected sequence: {}".format(len(sequence))
        return message
    else:
        # Return an empty string if data or value is None
        return ""


@app.callback(
    Output('seq-barplot', 'figure'),
    Input('dropdown', 'value')
)
def update_seq_barplot(value):
    # Check if value is None
    if value is None:
        # Return an empty string if value is None
        return ""

    # Calculate the frequency of each nucleotide in the selected sequence
    counts = {'A': value.count('A'), 'T': value.count('T'), 'G': value.count('G'),
              'C': value.count('C')}
    total_counts = sum(counts.values())
    try:
        freqs = {nt: count / total_counts for nt, count in counts.items()}
    except ZeroDivisionError:
        freqs = {nt: 0 for nt in counts}
    # Create a new bar plot with the updated data
    fig = go.Figure(data=[go.Bar(x=list(freqs.keys()), y=list(freqs.values()))],
                    layout=go.Layout(title='A T G C content', width=700, height=400))
    # fig.update_layout(template='seaborn')

    return fig


@app.callback(
    Output('high-GC-subsequence', 'children'),
    Input('dropdown', 'value'),
)
def get_highest_gc_subsequence(value):
    if value is None:
        return ""
    max_gc_content = -1
    max_gc_subseq = None
    for i in range(len(value) - 9):
        subseq = value[i:i + 10]
        gc_content = (subseq.count('G') + subseq.count('C')) / len(subseq)
        if gc_content > max_gc_content:
            max_gc_content = gc_content
            max_gc_subseq = subseq
    output = "Subsequence with high GC content: {}".format(max_gc_subseq)
    return output

@app.callback(
    Output('search-output', 'children'),
    Input('search-button', 'n_clicks'),
    State('search-input', 'value'),
    State('dropdown', 'value')
)
def search_callback(n_clicks, search_input, dropdown_value):
    if n_clicks:
        if dropdown_value is not None:
            if search_input in dropdown_value:
                return f"Subsequence '{search_input}' found in selected sequence!"
            else:
                return f"Subsequence '{search_input}' not found in selected sequence"
        else:
            return "Please select a sequence from the dropdown."
    else:
        return ""


if __name__ == '__main__':
    app.run_server(debug=True)
