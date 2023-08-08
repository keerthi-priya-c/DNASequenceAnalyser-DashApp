# DNASequenceAnalyser-DashApp

App Title: DNA Sequence Analyser https://app-vuwr.onrender.com/

The app is a DNA sequence analyser, built using Python and several libraries, including Dash, Plotly, and Pandas. The user interface is designed using Bootstrap. The app is a simple web application that allows the user to upload a CSV or a FASTA file (check the sample data folder for the required sequence file format). After the file is uploaded, the user can select a sequence from the file using a dropdown menu. The app displays the selected sequence and generates a bar plot showing the proportion of A, T, G, and C nucleotides in the selected sequence. The user can also search for a subsequence within the selected sequence. If the subsequence is found, the app highlights the subsequence in the displayed sequence.
## Installation and modules needed

•	import io\
•	import base64\
•	import dash\
•	import pandas as pd\
•	from Bio import SeqIO\
•	from dash import Input, Output, State, dcc, html\
•	import dash_bootstrap_components as dbc\
•	import plotly.graph_objs as go\


## Usage and Features
To run the application, navigate to the deployed app - https://app-vuwr.onrender.com/

To run the application using **Python**, execute:
app.py
Then, navigate to ` http://127.0.0.1:8050/ ` in your web browser.

Features of App:
1.	File Select or Drop file to Upload Area – csv or fasta file format input files are only allowed.
2.	The no. of Sequences in input file is displayed.
3.	Dropdown Menu to select the sequences present in input file.
4.	Length of selected sequence is displayed.
5.	The bar plot is populated reporting the proportion of A, T, G and C.
6.	Subsequence of length 10 with the highest G+C content is displayed.
7.	Additional Feature: Substring query- Any subsequence can be searched in selected sequence by entering the query and clicking search button.
 
