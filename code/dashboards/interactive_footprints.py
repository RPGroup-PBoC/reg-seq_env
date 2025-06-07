import numpy as np 
import pandas as pd 
import bokeh.io
import bokeh.plotting
from bokeh.models import * 
from bokeh.themes import Theme
from bokeh.transform import linear_cmap

import ast

import scipy

import glob
#import ast

#############################################
# Helper functions. Credit to Griffin Chure #
#############################################

def color_palette():
    """
    Returns a dictionary of the PBOC color palette
    """
    return {'green': '#7AA974', 'light_green': '#BFD598',
              'pale_green': '#DCECCB', 'yellow': '#EAC264',
              'light_yellow': '#F3DAA9', 'pale_yellow': '#FFEDCE',
              'blue': '#738FC1', 'light_blue': '#A9BFE3',
              'pale_blue': '#C9D7EE', 'red': '#D56C55', 'light_red': '#E8B19D',
              'pale_red': '#F1D4C9', 'purple': '#AB85AC',
              'light_purple': '#D4C2D9', 'dark_green':'#7E9D90', 'dark_brown':'#905426'}

def bokeh_theme():
    """A custom bokeh theme to match PBoC 2e colors"""
    theme_json = {'attrs':
            {'figure': {
                'background_fill_color': '#E3E7E9',
                'outline_line_color': '#FFFFFF',
            },
            'Axis': {
            'major_tick_in': 4,
            'major_tick_line_width': 1,
            'axis_label_text_font': 'Lato',
            'axis_label_text_font_style': 'normal'
            },
            'Grid': {
                'grid_line_color': "white",
            },
            'Legend': {
                'background_fill_color': '#E3E7E9',
                'border_line_color': '#FFFFFF',
                'border_line_width': 1.5,
                'background_fill_alpha': 0.5
            },
            'Text': {
                'text_font_style': 'normal',
               'text_font': 'Lato'
            },
            'Title': {
                'background_fill_color': '#FFFBCE',
                'text_font_style': 'normal',
                'align': 'center',
                'text_font': 'Lato',
                'offset': 2,
            }}}

    theme = Theme(json=theme_json)
    bokeh.io.curdoc().theme = theme

    # Define the colors
    colors = color_palette()
    palette = [v for k, v in colors.items() if 'pale' not in k]
    return [colors, palette]



def load_js(fname, args):
    """
    Given external javascript file names and arguments, load a bokeh CustomJS
    object
    
    Parameters
    ----------
    fname: str or list of str
        The file name of the external javascript file. If the desired javascript
        exists in multiple external files, they can be provided as a list of
        strings.
    args: dict
        The arguments to supply to the custom JS callback. 
    
    Returns
    -------
    cb : bokeh CustomJS model object
        Returns a bokeh CustomJS model object with the supplied code and
        arguments. This can be directly assigned as callback functions.
    """
    if type(fname) == str:
        with open('./js_scripts/'+fname) as f:
            js = f.read() 
    elif type(fname) == list:
        js = ''
        for _fname in fname:
            with open('./js_scripts/'+_fname) as f:
                js += f.read()

    cb = CustomJS(code=js, args=args)
    return cb

bokeh_theme()

# Path to store html file
bokeh.io.output_file('interactive_footprints.html')



def gaussian_smoothing(mut_info, sigma):
    return scipy.ndimage.gaussian_filter1d(mut_info, sigma)

def apply_window(pos, mut_info, smoothing_type, smoothing_kwargs):
    if smoothing_type == "gaussian":
        return pos, gaussian_smoothing(mut_info, smoothing_kwargs['sigma'])
    elif smoothing_type == "flat":
        d = smoothing_kwargs['d']
        filtered = np.array([np.mean(mut_info[i-d:i+d+1]) for i in np.arange(d, len(pos)-d)])
        return pos[d:-d], filtered

################################
# Data Import and Manipulation #
################################
collapse_df = pd.DataFrame()


df_gc = pd.read_csv("../../data/metadata/growth_conditions_short.csv", delimiter=';')

# Go through all files and compact the data
for file in glob.glob("../../data/footprints/*"):
    
    df = pd.read_csv(file)
    for promoter, group in df.groupby('promoter'):
        gc = file.split('/')[-1].split('-')[0]
        rep = file.split('/')[-1].split('_')[0][-1]
        x = group.pos.values
        y = group.mut_info.values
        collapse_df = pd.concat([collapse_df, pd.DataFrame(data={
            'mut_info':[y], 
            'pos': [x], 
            'promoter': promoter,
            'replicate': rep,
            'growth_condition': df_gc[df_gc['Index'] == int(gc)]['Condition'].values[0],
            })])
        
collapse_df = collapse_df.loc[[r not in ['ybeDp2', 'galEp'] for r in collapse_df.promoter], :]
#collapse_df = collapse_df[[x in ['araBp', 'araCp'] for x in collapse_df['promoter']]]
#collapse_df = collapse_df[[x in ['arabinose', 'glucose'] for x in collapse_df['growth_condition']]]

collapse_exshift_df = pd.DataFrame()


for file in glob.glob("../../data/expression_shifts/*"):
    df = pd.read_csv(file)
    for promoter, group in df.groupby('promoter'):
        pos = group.pos.values 
        base = group.base.values                                   
        wt_base = group.wt_base.values
        gc = file.split('/')[-1].split('-')[0]
        rep = file.split('/')[-1].split('_')[0][-1]
        expression_shift = group.expression_shift.values
        collapse_exshift_df = pd.concat([collapse_exshift_df, pd.DataFrame(data={
            'pos':[pos], 
            'base': [base], 
            'wt_base': [wt_base],
            'expression_shift': [expression_shift / (np.max(np.abs(expression_shift)))],
            'promoter': promoter,
            'replicate': rep,
            'growth_condition': df_gc[df_gc['Index'] == int(gc)]['Condition'].values[0]
            })])

collapse_exshift_df = collapse_exshift_df.loc[[r not in ['ybeDp2', 'galEp'] for r in collapse_exshift_df.promoter], :]
#collapse_exshift_df = collapse_exshift_df[[x in ['araBp', 'araCp'] for x in collapse_exshift_df['promoter']]]
#collapse_exshift_df = collapse_exshift_df[[x in ['arabinose', 'glucose'] for x in collapse_exshift_df['growth_condition']]]

# Import metadata for promoters
#df_meta = pd.read_csv('./20230525_footprints_meta.csv')
df_meta = pd.read_csv('./footprints_meta.csv')
df_regulonDB = pd.read_csv('./regulonDB_meta.csv')


def parse_space_separated_array(s):
    # Remove surrounding square brackets:
    s = s.strip().lstrip('[').rstrip(']')
    # Now parse the remaining space-separated values as floats:
    # (You can also use np.fromstring with sep=' ', but this split approach is simpler to debug.)
    return np.array([float(x) for x in s.split()])

df_hmm = pd.read_csv('./hmm_results.csv', converters={
                     "state_sequence": parse_space_separated_array,
                 })
df_hmm['growth_condition'] =  [df_gc[df_gc['Index'] == int(x)]['Condition'].values[0] for x in df_hmm['gc']]


# Transform types to strings
collapse_df['replicate'] = collapse_df['replicate'].astype(str)
collapse_exshift_df['replicate'] = collapse_exshift_df['replicate'].astype(str)
df_hmm['replicate'] = df_hmm['rep'].astype(str)

# Initiate CDS
data = ColumnDataSource(collapse_df)
exshift = ColumnDataSource(collapse_exshift_df)
meta = ColumnDataSource(df_meta)
regulonDB = ColumnDataSource(df_regulonDB)
promoters = list(df['promoter'].unique())
hmm = ColumnDataSource(df_hmm)


# Set inital settings for plot
prom_ini = 'araBp'
gc_ini = 'arabinose'
d_ini = 2
smoothing_type = "gaussian"
sigma_ini = 2



wt_seq = df_meta.loc[df_meta['promoter'] == prom_ini, "promoter_seq"].values[0]

# populate datasources with initial values
_df = collapse_df.loc[(collapse_df['promoter'] == prom_ini) 
           & (collapse_df['growth_condition'] == gc_ini),
           ['pos', 'mut_info', 'replicate']]


# put values in ColumnDataSource 
data_display = ColumnDataSource()

for rep, gdf in _df.groupby("replicate"):
    pos, mut_info = apply_window(
        pos=_df.loc[_df['replicate'] == rep, 'pos'].values[0], 
        mut_info=_df.loc[_df['replicate'] == rep, 'mut_info'].values[0],
        smoothing_type="gaussian",
        smoothing_kwargs={'sigma': sigma_ini}
    )
    
    data_display.data[f'pos_{rep}_1'] = pos
    data_display.data[f'mut_info_{rep}_1'] = mut_info
    a = np.zeros_like(pos, dtype=float)
    a[:] = np.nan
    data_display.data[f'pos_{rep}_2'] = a
    data_display.data[f'mut_info_{rep}_2'] =  a





# extract initial values for expression shift
exshift_display = ColumnDataSource()
_df_exshift = collapse_exshift_df.loc[(collapse_exshift_df['promoter'] == prom_ini) 
                                        & (collapse_exshift_df['growth_condition'] == gc_ini),
                                            ['pos', 'base', 'wt_base', 'expression_shift', 'replicate']]

for rep, gdf in _df_exshift.groupby("replicate"):
    exshift_display.data[f'pos_{rep}'] = gdf['pos'].values[0]
    exshift_display.data[f'base_{rep}'] = gdf['base'].values[0]
    exshift_display.data[f'wt_base_{rep}'] = gdf['wt_base'].values[0]
    exshift_display.data[f'expression_shift_{rep}'] = gdf['expression_shift'].values[0]


# Add empty arrays if 3rd replicate does not exist
if 'rep_3' not in data_display.data:
    data_display.data[f'pos_3_1'] = data_display.data['pos_1_1']
    data_display.data[f'pos_3_2'] = data_display.data['pos_1_2']
    data_display.data[f'mut_info_3_1'] = np.zeros_like(data_display.data['mut_info_1_1'])
    data_display.data[f'mut_info_3_2'] = np.zeros_like(data_display.data['mut_info_1_2'])
    
    exshift_display.data[f'pos_3'] = exshift_display.data['pos_1']
    exshift_display.data[f'base_3'] = exshift_display.data['base_1']
    exshift_display.data[f'wt_base_3'] = exshift_display.data['wt_base_1']
    exshift_display.data[f'expression_shift_3'] = np.zeros_like(exshift_display.data['expression_shift_1'])


# Compute cv for chosen data and fill CDS
cds_cv_gc = ColumnDataSource()
cds_cv_point_gc = ColumnDataSource()
cds_cv_patches = ColumnDataSource()

cds_cv_prom = ColumnDataSource()
cds_cv_point_prom = ColumnDataSource()

def cv(x):
    x = [gaussian_smoothing(y, sigma_ini) for y in x]
    return [np.std(y)/np.mean(y) for y in x]

# compute coefficient of variation for each footprint and store in CDS
df_cv_gc = collapse_df.loc[(collapse_df['promoter'] == prom_ini), :].groupby(["growth_condition", "replicate"])['mut_info'].agg(cv=cv).reset_index()
df_cv_gc.sort_values('cv', inplace=True)
df_cv_gc.reset_index(drop=True, inplace=True)
cds_cv_gc.data['x'] = [x[0] for x in df_cv_gc['cv'].values]
cds_cv_gc.data['y'] = np.arange(1, len(df_cv_gc)+1) / len(df_cv_gc)


cds_cv_patches.data['x_cv_1'] = [np.min(cds_cv_gc.data[f'x']), np.min(cds_cv_gc.data[f'x']), 0.65, 0.65]
cds_cv_patches.data['y'] = [0, 1, 1, 0]
cds_cv_patches.data['x_cv_2'] = [0.65, 0.65, 0.75, 0.75]
cds_cv_patches.data['x_cv_3'] = [0.75, 0.75, np.max(cds_cv_gc.data[f'x']), np.max(cds_cv_gc.data[f'x'])]

for rep in range(1, 4):
    indexes = df_cv_gc.loc[(df_cv_gc['growth_condition'] == gc_ini) & (df_cv_gc['replicate'] == f'{rep}'), :].index
    if len(indexes) == 0:
          cds_cv_point_gc.data[f'x_{rep}'] = [np.nan]
          cds_cv_point_gc.data[f'y_{rep}'] = [np.nan]
    else:
        ind = indexes[0]
        cds_cv_point_gc.data[f'x_{rep}'] = [cds_cv_gc.data[f'x'][ind]]
        cds_cv_point_gc.data[f'y_{rep}'] = [cds_cv_gc.data[f'y'][ind]]




df_cv_prom = collapse_df.loc[(collapse_df['growth_condition'] == gc_ini), :].groupby(["promoter", "replicate"])['mut_info'].agg(cv=cv).reset_index()
df_cv_prom.sort_values('cv', inplace=True)
df_cv_prom.reset_index(drop=True, inplace=True)
for rep, gdf in df_cv_prom.groupby('replicate'):
    cds_cv_prom.data[f'x_{rep}'] = [x[0] for x in gdf['cv'].values]
    cds_cv_prom.data[f'y_{rep}'] = np.arange(1, len(gdf)+1) / len(gdf)

    ind = np.where(gdf['promoter'] == prom_ini)[0][0]
    cds_cv_point_prom.data[f'x_{rep}'] = [cds_cv_prom.data[f'x_{rep}'][ind]]
    cds_cv_point_prom.data[f'y_{rep}'] = [cds_cv_prom.data[f'y_{rep}'][ind]]

    cds_cv_patches.data[f'x_prom_1_{rep}'] = [np.min(cds_cv_prom.data[f'x_{rep}']), np.min(cds_cv_prom.data[f'x_{rep}']), 0.65, 0.65]
    cds_cv_patches.data['y'] = [0, 1, 1, 0]
    cds_cv_patches.data[f'x_prom_2_{rep}'] = [0.65, 0.65, 0.75, 0.75]
    cds_cv_patches.data[f'x_prom_3_{rep}'] = [0.75, 0.75, np.max(cds_cv_prom.data[f'x_{rep}']), np.max(cds_cv_prom.data[f'x_{rep}'])]



# Fix third replicate   
if 'x_3' not in cds_cv_prom.data:

    a = np.empty(len(cds_cv_prom.data['x_1']))
    a[:] = np.nan
    cds_cv_prom.data['x_3'] = a
    cds_cv_prom.data['y_3'] = a
    cds_cv_point_prom.data['x_3'] = [np.nan]
    cds_cv_point_prom.data['y_3'] = [np.nan]
    cds_cv_patches.data['x_prom_1_3'] = [np.nan, np.nan, np.nan, np.nan]
    cds_cv_patches.data['x_prom_2_3'] = [np.nan, np.nan, np.nan, np.nan]
    cds_cv_patches.data['x_prom_3_3'] = [np.nan, np.nan, np.nan, np.nan]



###################
# Setting up plot #
###################

# Define the selections
prom_selector = Select(options=list(np.sort(promoters)), value=prom_ini)
gc_selector = Select(options=list(np.sort(collapse_df['growth_condition'].unique())), value=gc_ini)
sigma_slider = Slider(start=0.2, end=4, value=2, step=.2)
d_selector = Select(options=[str(x) for x in np.arange(6)], value=str(d_ini))
smooth_selector = Select(options=["gaussian", "flat"], value=str("gaussian"))
hmm_checkbox = Checkbox(label='Inferred Binding Sites', active=False)



# titles for selectors
prom_title = Div(text="<b>Promoter</b>")
gc_title = Div(text="<b>Growth Condition</b>")
smooth_title = Div(text="<b>Type of smoothing</b>")
sigma_title = Div(text="<b>Width of gaussian Kernel</b><br>(only if Gaussina smoothing is chosen)")
d_title = Div(text="<b>Window Width</b><br>(only if flat smoothing is chosen)")


# metadata for default choice
meta_ini = df_meta.loc[df_meta['promoter'] == prom_ini, :]


# boxes for description
prom_desc = Div(text='<div style="width:300px; overflow-wrap: break-word;"><b> Genes controlled by promoter</b>: <br/>' + meta_ini['genes'].values[0] + '<br/><b>Strand: </b><br/>' + meta_ini['direction'].values[0] + '<br/><b>5\':</b><br/>' + str(meta_ini['five_prime'].values[0]) + '<br/><b>3\':</b><br/>' + str(meta_ini['three_prime'].values[0]) + '</div>')
regulonDB_desc = Div(text="")

def update_sites(attr, old, new):
    x = '<div style="width:700px;"><b> Annotation in RegulonDB</b><br/>'
    promoter = prom_selector.value
    regulons = df_regulonDB.loc[df_regulonDB['PROMOTER_NAME'] == promoter, :]
    if len(regulons) == 0:
        x += '<br/>No Binding Sites Found'
    else:
        for index, site in regulons.iterrows():
            if -115 < site['CENTER_POSITION'] < 45:
                x += '<div style="overflow-wrap: break-word;"><br/><b>' + site['RI_FUNCTION'] + '</b><br/>Transcription Factor: ' + site['TRANSCRIPTION_FACTOR_NAME'] + '<br/>Binding Site Position Relative to TSS: ' + str(site['CENTER_POSITION']) + '</div>'#+ '<br/> Binding Site Sequence (Capital Letters): ' + site['RI_SEQUENCE'] + '<br/> Consensus Sequence: ' + site['CONSENSUS_SEQUENCE'] + '</div>';
    x += '</div>'
    regulonDB_desc.update(text=x)

update_sites("", "", "")




# Iterate through replicates and populate plot windows
TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

p_exshift = []
p_info = []
p_cvs = []

for i in range(1, 4):
    # expression shift plot
    p = bokeh.plotting.figure(width=1000, height=200, 
                                  x_axis_label='sequence',
                                  title="Expression Shift upon mutation",
                                  x_range=[-0.5-115, 45 - 0.5], 
                                  y_range=[0.5, 5 - 0.5],
                                  tooltips=[('wild type base', f'@wt_base_{i}')],
                                  tools=TOOLS)
    
    r = p.rect(x=f'pos_{i}', 
               y=f'base_{i}',
               width=1,
               height=1,
               fill_color=linear_cmap(f'expression_shift_{i}', 
                                      bokeh.palettes.interp_palette(["#63abbd", "#FFFFFF", "#ef865d"], 100),
                                      low=-1, 
                                      high=1
                                     ),
               line_color=None,
               source=exshift_display
    )   
    
    p.yaxis.ticker = np.arange(1,5)
    p.yaxis.major_label_overrides = {(tick+1): x_ for tick, x_ in enumerate(['A', 'C', 'G', 'T'])}

    p.xaxis.major_label_overrides = {(tick-115): x_ for tick, x_ in enumerate(exshift_display.data[f'wt_base_{i}'][0::4])}
    p.xaxis.major_label_text_font_size = "6pt"

    p.xaxis.ticker = np.arange(-115, 45)

    p.extra_x_ranges['x_above'] = Range1d(-115, 45)
    p.add_layout(LinearAxis(x_range_name='x_above', ticker=np.arange(-11, 5) * 10), 'above')

    color_bar = r.construct_color_bar(padding=5)
    p.add_layout(color_bar, "right")
    p_exshift.append(p)

    # footprint plot
    p = bokeh.plotting.figure(width=1000, height=200, 
                               x_axis_label='position',
                               y_axis_label='mutual information [bits]',
                               title="Mutual Information from Data")

    p.vbar(x=f'pos_{i}_1', top=f'mut_info_{i}_1', source=data_display)
    p.vbar(x=f'pos_{i}_2', top=f'mut_info_{i}_2', source=data_display, color='orange')

    p.xaxis.ticker = np.arange(-11, 5) * 10
    p_info.append(p)

    # ecdfs for cv
    p1 = bokeh.plotting.figure(width=300, height=300, 
                                  x_axis_label='coefficient of variation',
                                  y_axis_label='ECDF',
                                  title="Compared to all promoters in condition")
    p1.line(x=f'x_{i}', y=f'y_{i}', source=cds_cv_prom, line_width=2)
    p1.patch(x=f'x_prom_1_{i}', y='y', fill_color="#e59c8c", level='underlay', source=cds_cv_patches)
    p1.patch(x=f'x_prom_2_{i}', y='y', fill_color="#cdd6d1", level='underlay', source=cds_cv_patches)
    p1.patch(x=f'x_prom_3_{i}', y='y', fill_color="#bfd598", level='underlay', source=cds_cv_patches)
    p1.scatter(x=f'x_{i}', y=f'y_{i}', source=cds_cv_point_prom, size=5, line_width=1, fill_color="white", line_color="black")




    p2 = bokeh.plotting.figure(width=300, height=300, 
                                  x_axis_label='coefficient of variation',
                                  y_axis_label='ECDF',
                                  title="Compared to all conditions for promoter")
    p2.line(x=f'x', y=f'y', source=cds_cv_gc, line_width=2)
    p2.patch(x='x_cv_1', y='y', fill_color="#e59c8c", level='underlay', source=cds_cv_patches)
    p2.patch(x='x_cv_2', y='y', fill_color="#cdd6d1", level='underlay', source=cds_cv_patches)
    p2.patch(x='x_cv_3', y='y', fill_color="#bfd598", level='underlay', source=cds_cv_patches)
    p2.scatter(x=f'x_{i}', y=f'y_{i}', source=cds_cv_point_gc, size=5, line_width=1, fill_color="white", line_color="black")
    
    p_cvs.append([p1, p2])

# Define the callbacks
args = {
    'data_display': data_display,
    'exshift_display': exshift_display,
    'data': data,
    'exshift': exshift,
    'prom_selector': prom_selector,
    'gc_selector': gc_selector,
    'd_selector': d_selector,
    'sigma_slider': sigma_slider,
    'smooth_selector': smooth_selector,
    'prom_desc': prom_desc,
    'meta': meta,
    'regulonDB_desc': regulonDB_desc,
    'regulonDB': regulonDB,
    'x_axis': [p_exshift[i].xaxis[1] for i in range(3)],
    'p': p_exshift,
    'cds_cv_gc': cds_cv_gc,
    'cds_cv_prom': cds_cv_prom,
    'cds_cv_point_gc': cds_cv_point_gc,
    'cds_cv_point_prom': cds_cv_point_prom,
    'cds_cv_patches': cds_cv_patches,
    'p_cvs': p_cvs,
    'hmm_checkbox': hmm_checkbox,
    'hmm': hmm,
}




prom_cb = load_js(['data.js', 'functions.js', 'prom_selector.js', 'footprint_selector.js', 'cv.js'], args=args)
prom_selector.js_on_change('value', prom_cb)

gc_cb =  load_js(['data.js', 'functions.js', 'footprint_selector.js', 'cv.js'], args=args)
gc_selector.js_on_change('value', gc_cb)

footprint_cb = load_js(['data.js', 'functions.js',  'toggle_bs_off.js', 'footprint_selector.js','cv.js'], args=args)
for s in [smooth_selector, d_selector, sigma_slider, sigma_slider]:
    s.js_on_change('value', footprint_cb)

checkbox_cb =  load_js(['set_params_checkbox.js', 'data.js', 'functions.js', 'footprint_selector.js'], args=args)
hmm_checkbox.js_on_change('active', checkbox_cb)
selector_box = bokeh.layouts.row(
    bokeh.layouts.column(
        prom_title, 
        prom_selector,
        gc_title,
        gc_selector,
        hmm_checkbox
        
    ),
    bokeh.layouts.column( 
        smooth_title,
        smooth_selector,
        bokeh.layouts.row(
            bokeh.layouts.column(sigma_title,
                                 sigma_slider
            ),
            bokeh.layouts.column(d_title,
                                 d_selector
            )
        )
    ),
    prom_desc
)

plot_column = [bokeh.layouts.row(bokeh.layouts.column(p_info[i], p_exshift[i]), *(p_cvs[i])) for i in range(3)]

plot = bokeh.layouts.column(
    bokeh.layouts.column(
        selector_box,
        *plot_column
    ),
regulonDB_desc)

bokeh.io.save(plot)

# Remove first line from html document
with open(r'interactive_footprints.html', 'r+') as fp:
    # read an store all lines into list
    lines = fp.readlines()
    # move file pointer to the beginning of a file
    fp.seek(0)
    # truncate the file
    fp.truncate()

    # start writing lines except the first line
    # lines[1:] from line 2 to last line
    fp.writelines(lines[1:])


