import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
from bokeh.models import *
from bokeh.themes import Theme
from bokeh.transform import linear_cmap
import scipy.ndimage

#############################################
# Helper functions. Credit to Griffin Chure #
#############################################

def color_palette():
    return {'green': '#7AA974', 'light_green': '#BFD598',
              'pale_green': '#DCECCB', 'yellow': '#EAC264',
              'light_yellow': '#F3DAA9', 'pale_yellow': '#FFEDCE',
              'blue': '#738FC1', 'light_blue': '#A9BFE3',
              'pale_blue': '#C9D7EE', 'red': '#D56C55', 'light_red': '#E8B19D',
              'pale_red': '#F1D4C9', 'purple': '#AB85AC',
              'light_purple': '#D4C2D9', 'dark_green':'#7E9D90', 'dark_brown':'#905426'}

def bokeh_theme():
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
    colors = color_palette()
    palette = [v for k, v in colors.items() if 'pale' not in k]
    return [colors, palette]

bokeh_theme()
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

# Read the new single-file format
df_mi = pd.read_csv("mi_all.csv")
df_ex = pd.read_csv("ex_all.csv")

# Growth condition name lookup
df_gc = pd.read_csv("../../data/metadata/growth_conditions_short.csv", delimiter=';')

# Fix mislabeled conditions: swap 20-1 and 30-1
for df in [df_mi, df_ex]:
    mask_20 = df['ind'] == '20-1'
    mask_30 = df['ind'] == '30-1'
    df.loc[mask_20, 'ind'] = '__temp__'
    df.loc[mask_30, 'ind'] = '20-1'
    df.loc[df['ind'] == '__temp__', 'ind'] = '30-1'

# Parse ind column into gc and replicate
df_mi[['gc', 'replicate']] = df_mi['ind'].str.split('-', expand=True)
df_ex[['gc', 'replicate']] = df_ex['ind'].str.split('-', expand=True)

# Map gc index to condition name
gc_map = dict(zip(df_gc['Index'].astype(str), df_gc['Condition']))
df_mi['growth_condition'] = df_mi['gc'].map(gc_map)
df_ex['growth_condition'] = df_ex['gc'].map(gc_map)

# Shift MI positions from 0-159 to -115 to 44
df_mi['pos'] = df_mi['pos'] - 115

# Collapse MI data: group by promoter/gc/replicate, store arrays
collapse_df_fixed = []
for (prom, gc, rep), group in df_mi.groupby(['promoter', 'growth_condition', 'replicate']):
    group = group.sort_values('pos')
    collapse_df_fixed.append({
        'mut_info': group.mut_info.values * 1000,  # convert to mbits
        'pos': group.pos.values,
        'N': group.N.values,
        'promoter': prom,
        'replicate': rep,
        'growth_condition': gc
    })
collapse_df = pd.DataFrame(collapse_df_fixed)


# Collapse expression shift data
collapse_exshift_fixed = []
for (prom, gc, rep), group in df_ex.groupby(['promoter', 'growth_condition', 'replicate']):
    group = group.sort_values(['pos', 'base'])
    collapse_exshift_fixed.append({
        'pos': group.pos.values,
        'base': group.base.values,
        'wt_base': group.wt_base.values,
        'expression_shift': group.expression_shift.values,
        'promoter': prom,
        'replicate': rep,
        'growth_condition': gc
    })
collapse_exshift_df = pd.DataFrame(collapse_exshift_fixed)

# Compute average expression shift per position for coloring footprint bars
# Also normalize expression shifts per replicate for the heatmap
avg_ex_list = []
exshift_max_abs_list = []
for idx, row in collapse_exshift_df.iterrows():
    ex = row['expression_shift']
    max_abs = np.max(np.abs(ex))
    exshift_max_abs_list.append(max_abs if max_abs > 0 else 1.0)
    if max_abs > 0:
        collapse_exshift_df.at[idx, 'expression_shift'] = ex / max_abs
    # Compute mean expression shift per position (average over 4 bases)
    n_pos = len(ex) // 4
    avg_ex = np.array([np.mean(ex[i*4:(i+1)*4]) for i in range(n_pos)])
    avg_ex_list.append(avg_ex)
collapse_exshift_df['avg_expression_shift'] = avg_ex_list
collapse_exshift_df['max_abs_expression_shift'] = exshift_max_abs_list

# Import metadata
df_meta = pd.read_csv('./footprints_meta.csv')
df_regulonDB = pd.read_csv('./regulonDB_meta.csv')


# Initiate CDS
data = ColumnDataSource(collapse_df)
exshift = ColumnDataSource(collapse_exshift_df)
meta = ColumnDataSource(df_meta)
regulonDB = ColumnDataSource(df_regulonDB)
promoters = list(collapse_df['promoter'].unique())

# Set initial settings
prom_ini = promoters[0]
gc_options = list(np.sort(collapse_df.loc[collapse_df['promoter'] == prom_ini, 'growth_condition'].unique()))
gc_ini = gc_options[0]
d_ini = 2
smoothing_type = "gaussian"
sigma_ini = 2

# Initial data
_df = collapse_df.loc[(collapse_df['promoter'] == prom_ini) 
           & (collapse_df['growth_condition'] == gc_ini),
           ['pos', 'mut_info', 'replicate', 'N']]

available_reps = sorted(_df['replicate'].unique())
n_reps = len(available_reps)
has_third_rep = n_reps >= 3

# Build display CDS — slots are always 1, 2, 3 regardless of actual rep numbers
data_display = ColumnDataSource()

for slot, rep in enumerate(available_reps[:3], start=1):
    rep_data = _df.loc[_df['replicate'] == rep]
    pos_raw = rep_data['pos'].values[0]
    mi_raw = rep_data['mut_info'].values[0]
    N_raw = rep_data['N'].values[0]
    pos, mut_info = apply_window(pos_raw, mi_raw, "gaussian", {'sigma': sigma_ini})
    data_display.data[f'pos_{slot}'] = pos
    data_display.data[f'mut_info_{slot}'] = mut_info
    _, N_smooth = apply_window(pos_raw, N_raw, "gaussian", {'sigma': sigma_ini})
    data_display.data[f'N_{slot}'] = N_smooth
    data_display.data[f'color_{slot}'] = np.full(len(pos), '#738FC1')

# Fill remaining empty slots
for slot in range(n_reps + 1, 4):
    ref_pos = data_display.data['pos_1']
    data_display.data[f'pos_{slot}'] = ref_pos
    data_display.data[f'mut_info_{slot}'] = np.zeros(len(ref_pos))
    data_display.data[f'N_{slot}'] = np.ones(len(ref_pos))
    data_display.data[f'color_{slot}'] = np.full(len(ref_pos), '#738FC1')


# Expression shift display
exshift_display = ColumnDataSource()
_df_exshift = collapse_exshift_df.loc[(collapse_exshift_df['promoter'] == prom_ini) 
                                        & (collapse_exshift_df['growth_condition'] == gc_ini),
                                            ['pos', 'base', 'wt_base', 'expression_shift', 'replicate']]

available_reps_ex = sorted(_df_exshift['replicate'].unique())

for slot, rep in enumerate(available_reps_ex[:3], start=1):
    rep_data = _df_exshift.loc[_df_exshift['replicate'] == rep]
    exshift_display.data[f'pos_{slot}'] = rep_data['pos'].values[0]
    exshift_display.data[f'base_{slot}'] = rep_data['base'].values[0]
    exshift_display.data[f'wt_base_{slot}'] = rep_data['wt_base'].values[0]
    exshift_display.data[f'expression_shift_{slot}'] = rep_data['expression_shift'].values[0]

for slot in range(len(available_reps_ex) + 1, 4):
    ref = exshift_display.data.get('pos_1', np.arange(-115, 45).repeat(4))
    exshift_display.data[f'pos_{slot}'] = ref
    exshift_display.data[f'base_{slot}'] = exshift_display.data.get('base_1', np.tile([1,2,3,4], 160))
    exshift_display.data[f'wt_base_{slot}'] = np.full(len(ref), ' ')
    exshift_display.data[f'expression_shift_{slot}'] = np.zeros(len(ref))



###################
# Setting up plot #
###################

# Selectors
prom_selector = Select(options=list(np.sort(promoters)), value=prom_ini)
gc_selector = Select(options=gc_options, value=gc_ini)
sigma_slider = Slider(start=0.2, end=4, value=2, step=.2)
d_selector = Select(options=[str(x) for x in np.arange(6)], value=str(d_ini))
smooth_selector = Select(options=["gaussian", "flat"], value="gaussian")
bg_checkbox = Checkbox(label='Background Subtraction (8130/N)', active=False)
color_checkbox = Checkbox(label='Color by Expression Shift', active=False)

# Titles
prom_title = Div(text="<b>Promoter</b>")
gc_title = Div(text="<b>Growth Condition</b>")
smooth_title = Div(text="<b>Type of smoothing</b>")
sigma_title = Div(text="<b>Width of Gaussian Kernel</b><br>(only if Gaussian smoothing is chosen)")
d_title = Div(text="<b>Window Width</b><br>(only if flat smoothing is chosen)")

# Metadata display
meta_ini = df_meta.loc[df_meta['promoter'] == prom_ini, :]
if len(meta_ini) > 0:
    prom_desc = Div(text='<div style="width:300px; overflow-wrap: break-word;"><b>Genes controlled by promoter</b>: <br/>' + str(meta_ini['genes'].values[0]) + '<br/><b>Strand: </b><br/>' + str(meta_ini['direction'].values[0]) + '<br/><b>5\':</b><br/>' + str(meta_ini['five_prime'].values[0]) + '<br/><b>3\':</b><br/>' + str(meta_ini['three_prime'].values[0]) + '</div>')
else:
    prom_desc = Div(text='<div style="width:300px;">No metadata found</div>')
regulonDB_desc = Div(text="")


# Plots
TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

p_exshift = []
p_info = []

for i in range(1, 4):
    # Expression shift heatmap
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

    # Information footprint plot
    p = bokeh.plotting.figure(width=1000, height=200, 
                               x_axis_label='position',
                               y_axis_label='mutual information [mbits]',
                               title="Mutual Information from Data")
    p.vbar(x=f'pos_{i}', top=f'mut_info_{i}', source=data_display, fill_color={'field': f'color_{i}'}, line_color={'field': f'color_{i}'})
    p.xaxis.ticker = np.arange(-11, 5) * 10
    p_info.append(p)


# JavaScript functions
js_functions = """
function getAllIndexes(arr, val) {
    var indices = [], i = -1;
    while ((i = arr.indexOf(val, i+1)) != -1){
        indices.push(i);
    }
    return indices;
}

function gaussianFilter1D(data, sigma, options) {
    options = options || {};
    var mode = options.mode || 'reflect';
    var truncate = options.truncate || 4.0;
    var cval = options.cval || 0;
    if (sigma <= 0) return data.slice();
    var radius = Math.ceil(truncate * sigma);
    var kernelSize = 2 * radius + 1;
    var kernel = new Array(kernelSize);
    var sum = 0;
    var sigma2 = sigma * sigma;
    for (var i = -radius; i <= radius; i++) {
        var value = Math.exp(-(i * i) / (2 * sigma2));
        kernel[i + radius] = value;
        sum += value;
    }
    for (var i = 0; i < kernelSize; i++) kernel[i] /= sum;
    function reflectIndex(i, length) {
        while (i < 0 || i >= length) {
            if (i < 0) i = -i - 1;
            else if (i >= length) i = 2 * length - i - 1;
        }
        return i;
    }
    var result = new Array(data.length);
    for (var i = 0; i < data.length; i++) {
        var acc = 0;
        for (var j = -radius; j <= radius; j++) {
            var dataIndex = reflectIndex(i + j, data.length);
            acc += data[dataIndex] * kernel[j + radius];
        }
        result[i] = acc;
    }
    return result;
}

function flat_average(pos, mut_info, d) {
    if (d == 0) return [pos, mut_info];
    var smoothed = [];
    var new_pos = [];
    for (var i = d; i < (pos.length - d); i++) {
        new_pos.push(pos[i]);
        var s = 0;
        for (var k = i - d; k <= i + d; k++) s += mut_info[k];
        smoothed.push(s / (2 * d + 1));
    }
    return [new_pos, smoothed];
}

function smooth_data(pos, values, smooth_selector, d_selector, sigma_slider) {
    if (smooth_selector.value == "flat") {
        return flat_average(pos, values, Number(d_selector.value));
    } else {
        return [pos, gaussianFilter1D(values, Number(sigma_slider.value))];
    }
}

function exshift_to_color(avg_ex, n_pos) {
    // Custom colormap matching:
    // 0.0 -> #63ACBE, 0.1 -> #63ACBE, 0.5 -> #C8C8C8, 0.9 -> #EF875F, 1.0 -> #EF875F
    // val in [-1,1] maps to t in [0,1] via t = (val+1)/2
    var stops = [
        {t: 0.0, r: 0x63, g: 0xAC, b: 0xBE},
        {t: 0.1, r: 0x63, g: 0xAC, b: 0xBE},
        {t: 0.5, r: 0xC8, g: 0xC8, b: 0xC8},
        {t: 0.9, r: 0xEF, g: 0x87, b: 0x5F},
        {t: 1.0, r: 0xEF, g: 0x87, b: 0x5F}
    ];
    function interp_cmap(t) {
        if (t <= stops[0].t) return stops[0];
        if (t >= stops[stops.length-1].t) return stops[stops.length-1];
        for (var k = 0; k < stops.length - 1; k++) {
            if (t >= stops[k].t && t <= stops[k+1].t) {
                var frac = (t - stops[k].t) / (stops[k+1].t - stops[k].t);
                return {
                    r: Math.round(stops[k].r + (stops[k+1].r - stops[k].r) * frac),
                    g: Math.round(stops[k].g + (stops[k+1].g - stops[k].g) * frac),
                    b: Math.round(stops[k].b + (stops[k+1].b - stops[k].b) * frac)
                };
            }
        }
        return stops[stops.length-1];
    }
    var colors = new Array(n_pos);
    var max_abs = 0;
    for (var i = 0; i < avg_ex.length; i++) {
        if (Math.abs(avg_ex[i]) > max_abs) max_abs = Math.abs(avg_ex[i]);
    }
    if (max_abs == 0) max_abs = 1;
    for (var i = 0; i < n_pos; i++) {
        var val = (i < avg_ex.length) ? avg_ex[i] / max_abs : 0;  // -1 to 1
        var t = (val + 1) / 2;  // map to 0-1
        var c = interp_cmap(t);
        colors[i] = '#' + c.r.toString(16).padStart(2,'0') + c.g.toString(16).padStart(2,'0') + c.b.toString(16).padStart(2,'0');
    }
    return colors;
}
"""

# Main update callback JS
js_update = js_functions + """
var display = data_display.data;
var ex_display = exshift_display.data;
var source_data = data.data;
var exshift_data = exshift.data;

var prom = prom_selector.value;
var gc = gc_selector.value;
var bg_sub = bg_checkbox.active;
var color_by_ex = color_checkbox.active;

// Find matching indices in source data
var prom_inds = getAllIndexes(source_data['promoter'], prom);
var gc_inds = getAllIndexes(source_data['growth_condition'], gc);
var overlap_inds = prom_inds.filter(value => gc_inds.includes(value));

// Expression shift indices
var prom_inds_ex = getAllIndexes(exshift_data['promoter'], prom);
var gc_inds_ex = getAllIndexes(exshift_data['growth_condition'], gc);
var overlap_inds_ex = prom_inds_ex.filter(value => gc_inds_ex.includes(value));

// Collect available replicate numbers and sort them
var rep_set = new Set();
for (var k = 0; k < overlap_inds.length; k++) {
    rep_set.add(source_data['replicate'][overlap_inds[k]]);
}
var available_reps = Array.from(rep_set).sort();

var rep_set_ex = new Set();
for (var k = 0; k < overlap_inds_ex.length; k++) {
    rep_set_ex.add(exshift_data['replicate'][overlap_inds_ex[k]]);
}
var available_reps_ex = Array.from(rep_set_ex).sort();

// Assign available replicates to plot slots 1, 2, 3
for (var slot = 1; slot <= 3; slot++) {
    var rep_label = (slot <= available_reps.length) ? available_reps[slot - 1] : null;
    var rep_label_ex = (slot <= available_reps_ex.length) ? available_reps_ex[slot - 1] : null;

    if (rep_label !== null) {
        var rep_inds = getAllIndexes(source_data['replicate'], rep_label);
        var filtered = overlap_inds.filter(value => rep_inds.includes(value));
        var idx = filtered[0];
        
        var pos_raw = source_data['pos'][idx];
        var mi_raw = source_data['mut_info'][idx].slice();
        var N_raw = source_data['N'][idx];
        
        // Background subtraction BEFORE smoothing
        if (bg_sub) {
            for (var j = 0; j < mi_raw.length; j++) {
                mi_raw[j] = Math.max(0, mi_raw[j] - 8130 / N_raw[j]);
            }
        }
        
        var result = smooth_data(pos_raw, mi_raw, smooth_selector, d_selector, sigma_slider);
        var pos = result[0];
        var mut_info = result[1];
        
        var result_N = smooth_data(pos_raw, N_raw, smooth_selector, d_selector, sigma_slider);
        var N_smooth = result_N[1];
        
        display['pos_' + slot] = pos;
        display['mut_info_' + slot] = mut_info;
        display['N_' + slot] = N_smooth;
        
        // Color by expression shift
        if (color_by_ex && rep_label_ex !== null) {
            var rep_inds_ex2 = getAllIndexes(exshift_data['replicate'], rep_label_ex);
            var filtered_ex2 = overlap_inds_ex.filter(value => rep_inds_ex2.includes(value));
            if (filtered_ex2.length > 0) {
                var avg_ex = exshift_data['avg_expression_shift'][filtered_ex2[0]];
                display['color_' + slot] = exshift_to_color(avg_ex, pos.length);
            } else {
                var colors = new Array(pos.length);
                for (var j = 0; j < pos.length; j++) colors[j] = '#738FC1';
                display['color_' + slot] = colors;
            }
        } else {
            var colors = new Array(pos.length);
            for (var j = 0; j < pos.length; j++) colors[j] = '#738FC1';
            display['color_' + slot] = colors;
        }
    } else {
        // Empty slot
        var ref_len = display['pos_1'].length;
        display['mut_info_' + slot] = new Array(ref_len).fill(0);
        var colors = new Array(ref_len);
        for (var j = 0; j < ref_len; j++) colors[j] = '#738FC1';
        display['color_' + slot] = colors;
    }
    
    // Expression shift display
    if (rep_label_ex !== null) {
        var rep_inds_ex3 = getAllIndexes(exshift_data['replicate'], rep_label_ex);
        var filtered_ex3 = overlap_inds_ex.filter(value => rep_inds_ex3.includes(value));
        if (filtered_ex3.length > 0) {
            var idx_ex = filtered_ex3[0];
            ex_display['pos_' + slot] = exshift_data['pos'][idx_ex];
            ex_display['base_' + slot] = exshift_data['base'][idx_ex];
            ex_display['wt_base_' + slot] = exshift_data['wt_base'][idx_ex];
            ex_display['expression_shift_' + slot] = exshift_data['expression_shift'][idx_ex];
        } else {
            ex_display['expression_shift_' + slot] = new Array(ex_display['expression_shift_' + slot].length).fill(0);
            ex_display['wt_base_' + slot] = new Array(ex_display['wt_base_' + slot].length).fill(" ");
        }
    } else {
        ex_display['expression_shift_' + slot] = new Array(ex_display['expression_shift_' + slot].length).fill(0);
        ex_display['wt_base_' + slot] = new Array(ex_display['wt_base_' + slot].length).fill(" ");
    }
    
    // Update x-axis labels
    var map1 = new Map();
    for (var ii = 0; ii < ex_display['wt_base_' + slot].length; ii += 4) {
        var j = (ii) / 4 - 115;
        map1.set(j, ex_display['wt_base_' + slot][ii]);
    }
    x_axis[slot - 1].major_label_overrides = map1;
    p[slot - 1].reset.emit();
    p[slot - 1].change.emit();
}

// Toggle third replicate row visibility
row_3.visible = (available_reps.length >= 3);

// Update plot titles with actual replicate numbers
for (var slot = 1; slot <= 3; slot++) {
    var rep_label = (slot <= available_reps.length) ? available_reps[slot - 1] : "";
    if (rep_label) {
        p_info[slot - 1].title.text = "Mutual Information (Replicate " + rep_label + ")";
        p[slot - 1].title.text = "Expression Shift (Replicate " + rep_label + ")";
    }
}

data_display.change.emit();
exshift_display.change.emit();
"""

# Promoter selector JS (updates gc_selector options + triggers main update)
js_prom_selector = js_functions + """
var prom = prom_selector.value;
var source_data = data.data;
var prom_inds = getAllIndexes(source_data['promoter'], prom);
var gc_set = new Set();
for (var i = 0; i < prom_inds.length; i++) {
    gc_set.add(source_data['growth_condition'][prom_inds[i]]);
}
var gc_list = Array.from(gc_set).sort();
gc_selector.options = gc_list;
if (!gc_list.includes(gc_selector.value)) {
    gc_selector.value = gc_list[0];
}

// Update metadata
var meta_inds = getAllIndexes(meta.data['promoter'], prom);
if (meta_inds.length > 0) {
    var direction = meta.data['direction'][meta_inds[0]];
    var genes = meta.data['genes'][meta_inds[0]];
    var five_prime = meta.data['five_prime'][meta_inds[0]];
    var three_prime = meta.data['three_prime'][meta_inds[0]];
    prom_desc.text = '<div style="width:300px; overflow-wrap: break-word;"><b>Genes controlled by promoter</b>: <br/>' + genes + '<br/><b>Strand: </b><br/>' + direction + '<br/><b>5\\':</b><br/>' + five_prime + '<br/><b>3\\':</b><br/>' + three_prime + '</div>';
}

// Update RegulonDB annotation
regulonDB_desc.text = '<div style="width:700px;"><b>Annotation in RegulonDB</b><br/>';
var regulon_indices = getAllIndexes(regulonDB.data['PROMOTER_NAME'], prom);
if (regulon_indices.length == 0) {
    regulonDB_desc.text += '<br/>No Binding Sites Found';
} else {
    for (var i = 0; i < regulon_indices.length; i++) {
        regulonDB_desc.text += '<div style="overflow-wrap: break-word;"><br/><b>' + regulonDB.data['RI_FUNCTION'][regulon_indices[i]] + '</b><br/>Transcription Factor: ' + regulonDB.data['TRANSCRIPTION_FACTOR_NAME'][regulon_indices[i]] + '<br/>Binding Site Position Relative to TSS: ' + regulonDB.data['CENTER_POSITION'][regulon_indices[i]] + '</div>';
    }
}
regulonDB_desc.text += '</div>';
"""

# Build the row for the 3rd replicate (will be toggled)
row_3 = bokeh.layouts.column(p_info[2], p_exshift[2], visible=has_third_rep)

# Download controls
download_data_selector = Select(
    options=["Mutual Information", "Expression Shifts", "Both"], 
    value="Mutual Information"
)
download_scope_selector = Select(
    options=["Current view", "All conditions for promoter", "All promoters for condition", "Everything"], 
    value="Current view"
)
download_button = Button(label="Download CSV", button_type="success")

download_title = Div(text="<b>Download</b>")
download_data_title = Div(text="Data type")
download_scope_title = Div(text="Scope")

# Download JS callback
js_download = js_functions + """
var source_data = data.data;
var exshift_data = exshift.data;
var prom = prom_selector.value;
var gc = gc_selector.value;
var data_type = download_data_selector.value;
var scope = download_scope_selector.value;

// Determine which indices to include
function get_indices(scope, prom, gc) {
    var prom_inds = getAllIndexes(source_data['promoter'], prom);
    var gc_inds = getAllIndexes(source_data['growth_condition'], gc);
    
    if (scope == "Current view") {
        return prom_inds.filter(v => gc_inds.includes(v));
    } else if (scope == "All conditions for promoter") {
        return prom_inds;
    } else if (scope == "All promoters for condition") {
        return gc_inds;
    } else {
        // Everything
        var all = [];
        for (var i = 0; i < source_data['promoter'].length; i++) all.push(i);
        return all;
    }
}

function get_indices_ex(scope, prom, gc) {
    var prom_inds = getAllIndexes(exshift_data['promoter'], prom);
    var gc_inds = getAllIndexes(exshift_data['growth_condition'], gc);
    
    if (scope == "Current view") {
        return prom_inds.filter(v => gc_inds.includes(v));
    } else if (scope == "All conditions for promoter") {
        return prom_inds;
    } else if (scope == "All promoters for condition") {
        return gc_inds;
    } else {
        var all = [];
        for (var i = 0; i < exshift_data['promoter'].length; i++) all.push(i);
        return all;
    }
}

var csv_parts = [];

function csvField(val) {
    var s = "" + val;
    if (s.indexOf(",") >= 0 || s.indexOf('"') >= 0) {
        return '"' + s.replace(/"/g, '""') + '"';
    }
    return s;
}

// Mutual information CSV
if (data_type == "Mutual Information" || data_type == "Both") {
    var mi_rows = ["promoter,growth_condition,replicate,pos,mut_info,N"];
    var inds = get_indices(scope, prom, gc);
    for (var k = 0; k < inds.length; k++) {
        var idx = inds[k];
        var p = source_data['promoter'][idx];
        var g = source_data['growth_condition'][idx];
        var r = source_data['replicate'][idx];
        var pos_arr = source_data['pos'][idx];
        var mi_arr = source_data['mut_info'][idx];
        var N_arr = source_data['N'][idx];
        for (var j = 0; j < pos_arr.length; j++) {
            mi_rows.push(csvField(p) + "," + csvField(g) + "," + r + "," + pos_arr[j] + "," + mi_arr[j] + "," + N_arr[j]);
        }
    }
    csv_parts.push({name: "mutual_information.csv", content: mi_rows.join("\\n")});
}

// Expression shift CSV
if (data_type == "Expression Shifts" || data_type == "Both") {
    var ex_rows = ["promoter,growth_condition,replicate,pos,base,wt_base,expression_shift"];
    var inds_ex = get_indices_ex(scope, prom, gc);
    for (var k = 0; k < inds_ex.length; k++) {
        var idx = inds_ex[k];
        var p = exshift_data['promoter'][idx];
        var g = exshift_data['growth_condition'][idx];
        var r = exshift_data['replicate'][idx];
        var pos_arr = exshift_data['pos'][idx];
        var base_arr = exshift_data['base'][idx];
        var wt_arr = exshift_data['wt_base'][idx];
        var ex_arr = exshift_data['expression_shift'][idx];
        var max_abs = exshift_data['max_abs_expression_shift'][idx];
        for (var j = 0; j < pos_arr.length; j++) {
            ex_rows.push(csvField(p) + "," + csvField(g) + "," + r + "," + pos_arr[j] + "," + base_arr[j] + "," + csvField(wt_arr[j]) + "," + (ex_arr[j] * max_abs));
        }
    }
    csv_parts.push({name: "expression_shifts.csv", content: ex_rows.join("\\n")});
}

// Trigger download(s)
for (var f = 0; f < csv_parts.length; f++) {
    var blob = new Blob([csv_parts[f].content], {type: 'text/csv'});
    var url = URL.createObjectURL(blob);
    var a = document.createElement('a');
    a.href = url;
    a.download = csv_parts[f].name;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}
"""

# Args for callbacks
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
    'bg_checkbox': bg_checkbox,
    'color_checkbox': color_checkbox,
    'prom_desc': prom_desc,
    'meta': meta,
    'regulonDB_desc': regulonDB_desc,
    'regulonDB': regulonDB,
    'x_axis': [p_exshift[i].xaxis[1] for i in range(3)],
    'p': p_exshift,
    'p_info': p_info,
    'row_3': row_3,
}

# Callbacks
update_cb = CustomJS(code=js_update, args=args)
prom_cb = CustomJS(code=js_prom_selector + "\n" + js_update, args=args)

prom_selector.js_on_change('value', prom_cb)
gc_selector.js_on_change('value', update_cb)
for s in [smooth_selector, d_selector, sigma_slider]:
    s.js_on_change('value', update_cb)
bg_checkbox.js_on_change('active', update_cb)
color_checkbox.js_on_change('active', update_cb)

download_args = {
    'data': data,
    'exshift': exshift,
    'prom_selector': prom_selector,
    'gc_selector': gc_selector,
    'download_data_selector': download_data_selector,
    'download_scope_selector': download_scope_selector,
}
download_cb = CustomJS(code=js_download, args=download_args)
download_button.js_on_click(download_cb)


# Layout
selector_box = bokeh.layouts.row(
    bokeh.layouts.column(
        prom_title, 
        prom_selector,
        gc_title,
        gc_selector,
        bg_checkbox,
        color_checkbox,
    ),
    bokeh.layouts.column( 
        smooth_title,
        smooth_selector,
        bokeh.layouts.row(
            bokeh.layouts.column(sigma_title, sigma_slider),
            bokeh.layouts.column(d_title, d_selector)
        )
    ),
    prom_desc,
    bokeh.layouts.column(
        download_title,
        download_data_title,
        download_data_selector,
        download_scope_title,
        download_scope_selector,
        download_button,
    ),
)

plot_rows = [
    bokeh.layouts.column(p_info[0], p_exshift[0]),
    bokeh.layouts.column(p_info[1], p_exshift[1]),
    row_3,
]

plot = bokeh.layouts.column(
    selector_box,
    *plot_rows,
    regulonDB_desc
)

bokeh.io.save(plot)

# Remove first line from html document
with open(r'interactive_footprints.html', 'r+') as fp:
    lines = fp.readlines()
    fp.seek(0)
    fp.truncate()
    fp.writelines(lines[1:])