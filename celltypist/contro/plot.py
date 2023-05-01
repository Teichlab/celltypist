import numpy as np
import pandas as pd
from typing import Union, Optional
from .symbols import NOVEL, REMAIN
import sys
try:
    import plotly.graph_objects as go
except ImportError:
    pass
import matplotlib
from matplotlib import pyplot as plt
import os
SEP1 = '_____'
SEP2 = 'CXCXCXCX'
DEFAULT_SANKEY_COLORS = ['#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#B68100', '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16',
                         '#A777F1', '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038', '#FD3216', '#00FE35', '#6A76FC', '#FED4C4', '#FE00CE', '#0DF9FF', '#F6F926', '#FF9616', '#479B55',
                         '#EEA6FB', '#DC587D', '#D626FF', '#6E899C', '#00B5F7', '#B68E00', '#C9FBE5', '#FF0092', '#22FFA7', '#E3EE9E', '#86CE00', '#BC7196', '#7E7DCD', '#FC6955', '#E48F72']

def _qq_order(sub_relation, ordered_cols):
    """
    Disentangle the cross connections in tree plot.
    """
    sub_relation2 = sub_relation[ordered_cols].copy()
    nms = sub_relation2.iloc[:, 0].value_counts().index.values
    for i in range(len(ordered_cols) - 1):
        i += 1
        val1 = sub_relation2.iloc[:, :i].agg('_'.join, axis = 1)
        val2 = sub_relation2.iloc[:, i]
        tb = pd.crosstab(val1, val2).loc[nms]
        nms = []
        flag = 0
        for j in range(tb.shape[0]):
            ind1 = (tb.iloc[j] > 0).values
            if j < (tb.shape[0] - 1):
                ind2 = (tb.iloc[j + 1] > 0).values
            else:
                ind2 = np.repeat([False], tb.shape[1])
            ind = ind1 & ind2
            all_nms = tb.columns[ind1][np.argsort(tb.loc[:, ind1].iloc[j])[::-1]]
            nm1 = tb.index[j]
            if ind.sum() > 0:
                last_nm = tb.columns[ind][0]
                mid_nms = np.setdiff1d(all_nms, last_nm)
                if flag == 1:
                    first_nm = fnm
                    mid_nms = np.setdiff1d(mid_nms, first_nm)
                    nms.append(nm1 + '_' + first_nm)
                for md in mid_nms:
                    nms.append(nm1 + '_' + md)
                nms.append(nm1 + '_' + last_nm)
                flag = 1
                fnm = last_nm
            else:
                if flag == 1:
                    all_nms = np.setdiff1d(all_nms, fnm)
                    nms.append(nm1 + '_' + fnm)
                for md in all_nms:
                    nms.append(nm1 + '_' + md)
                flag = 0
    sub_relation.index = sub_relation2.agg('_'.join, axis = 1)
    return sub_relation.loc[nms]

#def _qq_order2(sub_relation, ordered_cols):
#    sub_relation2 = sub_relation[ordered_cols].copy()
#    nms = sub_relation2.iloc[:, 0].value_counts().index.values
#    for i in range(len(ordered_cols)-1):
#        i += 1
#        val1 = sub_relation2.iloc[:, :i].agg('@'.join, axis=1)
#        val2 = sub_relation2.iloc[:, i]
#        tb = pd.crosstab(val1, val2)
#        tb = tb.loc[nms]
#        csum = (tb > 0).sum(0)
#        if csum.max() > 1:
#            idxs = np.argsort(csum)
#            for idx in idxs[csum[idxs] > 1]:
#                idxs2 = np.where(tb.iloc[:, idx] > 0)[0]
#                if np.any(np.diff(idxs2) > 1):
#                    ti = idxs2[0]
#                    tnm = np.array(nms[ti].split('@'))
#                    for idx2 in idxs2[1:]:
#                        nm = np.array(nms[idx2].split('@'))
#                        if idx2 == 0:
#                            nm2 = nms[idx2+1].split('_')
#                        elif idx2 + 1 < len(nms):
#                            nm2 = nms[idx2-1].split('_') + nms[idx2+1].split('_')
#                        else:
#                            nm2 = nms[idx2-1].split('_')
#                        if len(set(nm) & set(nm2)) == 0:
#                            mv_idxs.append(idx2)
#        nms = []
#        flag = 0
#        for j in range(tb.shape[0]):
#            ind1 = (tb.iloc[j] > 0).values
#            if j < (tb.shape[0] - 1):
#                ind2 = (tb.iloc[j+1] > 0).values
#            else:
#                ind2 = np.repeat([False], tb.shape[1])
#            ind = (ind1 & ind2)
#            all_nms = tb.columns[ind1][np.argsort(tb.loc[:, ind1].iloc[j])[::-1]]
#            nm1 = tb.index[j]
#            if ind.sum() > 0:
#                last_nm = tb.columns[ind][0]
#                mid_nms = np.setdiff1d(all_nms, last_nm)
#                if flag == 1:
#                    first_nm = fnm
#                    mid_nms = np.setdiff1d(mid_nms, first_nm)
#                    nms.append(nm1 + '_' + first_nm)
#                for md in mid_nms:
#                    nms.append(nm1 + '_' + md)
#                nms.append(nm1 + '_' + last_nm)
#                flag = 1
#                fnm = last_nm
#            else:
#                if flag == 1:
#                    all_nms = np.setdiff1d(all_nms, fnm)
#                    nms.append(nm1 + '_' + fnm)
#                for md in all_nms:
#                    nms.append(nm1 + '_' + md)
#                flag = 0
#    sub_relation.index = sub_relation2.agg('_'.join, axis=1)
#    sub_relation = sub_relation.loc[nms]
#    return sub_relation

def _relation_to_data(relation, return_sankey: bool = True) -> Union[pd.DataFrame, tuple]:
    """
    For internal use. Turn the harmonization result into Sankey input.
    """
    datasets = relation.columns[0::2]
    relation = relation[datasets].copy()
    relation.index = np.arange(relation.shape[0])
    #add prefix for cell types
    for i in relation.index:
        for j in datasets:
            content = relation.loc[i, j]
            if content not in [NOVEL, REMAIN]:
                relation.loc[i, j] = j + SEP1 + content
    #rename columns
    for _j in range(len(datasets)):
        j = datasets[_j]
        mapping_NOVEL = dict()
        mapping_REMAIN = dict()
        suffix_NOVEL = 0
        suffix_REMAIN = 0
        for i in relation.index:
            content = relation.loc[i, j]
            if content in [NOVEL, REMAIN]:
                celltypes = relation.loc[i].values
                relay = celltypes[~np.isin(celltypes, [NOVEL, REMAIN])][0] if _j == 0 else celltypes[_j-1]
                if content == NOVEL:
                    if relay not in mapping_NOVEL:
                        suffix_NOVEL += 1
                        mapping_NOVEL[relay] = f"{NOVEL}{SEP1}{suffix_NOVEL}"
                    relation.loc[i, j] = j + SEP1 + mapping_NOVEL[relay]
                else:
                    if relay not in mapping_REMAIN:
                        suffix_REMAIN += 1
                        mapping_REMAIN[relay] = f"{REMAIN}{SEP1}{suffix_REMAIN}"
                    relation.loc[i, j] = j + SEP1 + mapping_REMAIN[relay]
    if not return_sankey:
        return relation
    #label
    label = np.concatenate([np.unique(relation[j]) for j in datasets])
    refer = pd.Series(label).reset_index().set_index(0)
    #link = source + target + value
    link = pd.DataFrame(columns = ['source', 'target', 'value'])
    for k in range(len(datasets) - 1):
        vc = (relation[datasets[k]] + SEP2 + relation[datasets[k+1]]).value_counts().reset_index()
        vc[['source', 'target']] = vc['index'].str.split(SEP2, expand=True).values
        vc.rename(columns = {0: 'value'}, inplace = True)
        link = pd.concat([link, vc[['source', 'target', 'value']]])
    source = refer.loc[link.source.values, 'index'].values
    target = refer.loc[link.target.values, 'index'].values
    value = link.value.astype(int).values
    return relation, go.Sankey(node = dict(label = label), link = dict(source = source, target = target, value = value))

def _identify_relation_groups(relation, group_prefix: str = 'Group', order_row: bool = True, order_column: bool = False) -> tuple:
    """
    For internal use. Identify cell type groups based on the cell type harmonization result.
    """
    new_relation = _relation_to_data(relation, False)
    datasets = new_relation.columns
    dup_celltypes = np.unique(np.concatenate([new_relation[dataset].values[new_relation[dataset].duplicated()] for dataset in datasets]))
    if len(dup_celltypes) > 0:
        groups = np.full(new_relation.shape[0], f"{group_prefix}0", dtype = object)
        rownames = new_relation.index
        i = 0
        while len(rownames) > 0:
            i += 1
            receive = []
            provide = [rownames[0]]
            while len(provide) > 0:
                pp = provide.pop()
                receive.append(pp)
                for dataset in datasets:
                    celltype = new_relation.loc[pp, dataset]
                    if celltype in dup_celltypes:
                        extends = new_relation[new_relation[dataset] == celltype].index
                        provide.extend(extends.tolist())
                provide = list(set(provide).difference(receive))
            groups[new_relation.index.isin(receive)] = f"{group_prefix}{i}"
            rownames = rownames[~np.isin(rownames, receive)]
    else:
        groups = np.array([f"{group_prefix}{i+1}" for i in range(new_relation.shape[0])], dtype = object)
    ##remove in the future-->
    assert np.all(groups != f"{group_prefix}0")
    ##<<-remove in the future
    if order_row:
        df = []
        gs = []
        for j in range(1, len(np.unique(groups))+1):
            each_group = f"{group_prefix}{j}"
            sub_relation = new_relation[groups == each_group]
            if sub_relation.shape[0] > 1:
                col_uniques = sub_relation.apply(pd.Series.nunique).values
                ordered_cols = datasets[np.argsort(col_uniques)]
                #sub_relation = sub_relation.sort_values(by = ordered_cols.tolist())
                sub_relation = _qq_order(sub_relation, ordered_cols)
                if order_column:
                    sub_relation = sub_relation[ordered_cols]
                    sub_relation.columns = [f"D{x+1}" for x in range(sub_relation.shape[1])]
                    flag_unique = ~sub_relation.isin(dup_celltypes)
                    for row_index in range(flag_unique.shape[0]):
                        row_series = flag_unique.iloc[row_index]
                        col_indices = []
                        last_col = flag_unique.shape[1] - 1
                        while row_series.values[last_col]:
                            col_indices.append(last_col)
                            last_col -= 1
                        col_indices.reverse()
                        if len(col_indices) >= 2:
                            is_blank = np.isin([x.split(SEP1)[1] for x in sub_relation.iloc[row_index, col_indices]], [NOVEL, REMAIN])
                            sub_relation.iloc[row_index, col_indices] = sub_relation.iloc[row_index, col_indices].values[np.argsort(~is_blank)]
            else:
                if order_column:
                    is_blank = np.isin(relation.loc[groups == each_group, datasets].values[0], [NOVEL, REMAIN])
                    sub_relation = sub_relation[datasets[np.argsort(~is_blank)]]
                    sub_relation.columns = [f"D{x+1}" for x in range(sub_relation.shape[1])]
            df.append(sub_relation)
            gs.extend([each_group] * sub_relation.shape[0])
        return np.array(gs, dtype = object), pd.concat(df, axis = 0, ignore_index = True)
    return groups, new_relation

def _mix_colors(cols) -> str:
    """
    For internal use. Get the blended color.
    """
    return matplotlib.colors.to_hex(np.array([matplotlib.colors.to_rgb(col) for col in cols]).mean(axis = 0))

def _new_relation_to_color(new_relation, node_color, novel_node_color, remain_node_color, cmap = 'Reds') -> dict:
    """
    For internal use. Get the cell-type-to-color mapping.
    """
    map_color = {NOVEL: novel_node_color, REMAIN: remain_node_color}
    if node_color is None:
        for i in new_relation.index:
            row_color = DEFAULT_SANKEY_COLORS[i % len(DEFAULT_SANKEY_COLORS)]
            for j in new_relation.columns:
                content = new_relation.loc[i, j]
                if content.split(SEP1)[1] in [NOVEL, REMAIN]:
                    continue
                if content not in map_color:
                    map_color[content] = row_color
                else:
                    map_color[content] = _mix_colors([map_color[content], row_color])
    elif isinstance(node_color, pd.DataFrame):
        node_color = node_color.copy()
        map_values = node_color.iloc[:, 2].values
        if not isinstance(map_values[0], str):
            q10 = np.quantile(map_values, 0.10)
            q85 = np.quantile(map_values, 0.85)
            map_values[map_values <= q10] = q10
            map_values[map_values >= q85] = q85
            map_values = (map_values - map_values.min()) / map_values.ptp()
            map_values = plt.get_cmap(cmap, 256)(map_values)
            map_values = [matplotlib.colors.to_hex(map_value) for map_value in map_values]
        node_color['_combination'] = node_color.iloc[:, 0].astype(str) + SEP1 + node_color.iloc[:, 1].astype(str)
        map_color.update(dict(zip(node_color['_combination'].values, map_values)))
    else:
        raise TypeError(
                f"🛑 Please provide `node_color` as a data frame")
    return map_color

def sankeyplot(alignment,
           #node colors
           node_color: Optional[pd.DataFrame] = None, novel_node_color: str = '#FFFFFF', remain_node_color: str = '#F0F0F0',
           #link color
           link_color: Optional[str] = None,
           #figure elements
           title: str = 'CellTypist label harmonization',
           #figure size
           show: bool = True, save: Union[str, bool] = False, width: Optional[int] = None, height: Optional[int] = None,
           #to fig.update_layout
           layout_dict: dict = {},
           #to fig.update_traces
           trace_dict: dict = {},
           #for developer use
           expand_label: bool = False,
          ) -> None:
    """
    Generate a Sankey diagram showing the CellTypist label harmonization in a qualitative manner.

    Parameters
    ----------
    alignment
        A :class:`~celltypist.contro.align.DistanceAlignment` or :class:`~pandas.DataFrame` object representing the harmonization result.
    node_color
        A :class:`~pandas.DataFrame` with three consecutive columns representing dataset, cell type, and color, respectively.
        Default to a color scheme that allows matched cell types to have the same colors.
    novel_node_color
        Color of dataset-specific (i.e., novel) cell types.
        (Default: `'#FFFFFF'`)
    remain_node_color
        Color of remaining unresolved cell types.
        (Default: `'#F0F0F0'`)
    link_color
        Color of links. Default to translucent grey as used in plotly.
    title
        Figure title.
        (Default: `'CellTypist label harmonization'`)
    show
        Whether to show the plot.
        (Default: `True`)
    save
        Whether to save the plot. This can also be a figure filename.
        Supported figure suffixes are: .html, .png, .jpg, .jpeg, .webp, .svg, .pdf, .eps.
        (Default: `False`)
    width
        Figure width in pixels.
        Default to 700 in a canonical Plotly setting.
    height
        Figure height in pixels.
        Default to 450 in a canonical Plotly setting.
    layout_dict
        A dict passed to the method `.update_layout` of :class:`plotly.graph_objects.Figure` for setting the figure layout.
        Example keys include `plot_bgcolor` which sets the plot area color, `font_color` which sets the text color, `font_size` which sets the text size, etc.
    trace_dict
        A dict passed to the method `.update_traces` of :class:`plotly.graph_objects.Figure` for setting the Sankey plot.
        Example keys include `note_pad` which sets the padding between nodes, `link_line_color` which sets the link border color, `orientation` which sets the plot orientation, etc.
    expand_label
        Ignored. Whether to show the unique expanded labels. Only for developer use.
        (Default: `False`)

    Returns
    ----------
    None
    """
    if 'plotly' not in sys.modules:
        logger.warn(f"⚠️ Warning: to draw a Sankey diagram, package `plotly` is required. Please install `plotly` first")
        return
    if isinstance(alignment, pd.DataFrame):
        relation = alignment
    elif hasattr(alignment, 'relation'):
        relation = alignment.relation
    else:
        raise TypeError(
                f"🛑 Please provide correct input - either a DistanceAlignment or a data frame")
    trace = _relation_to_data(relation, return_sankey = True)[1]
    #relation2new_relation is run twice actually, but time cost is negligible; this new relation is row ordered
    new_relation = _identify_relation_groups(relation, group_prefix = 'Group', order_row = True, order_column = False)[1]
    expanded_label = trace.node.label
    original_label = pd.Series(expanded_label).str.split(SEP1, expand = True)[1].values
    blank_flag = np.isin(original_label, [NOVEL, REMAIN])
    #node color
    if isinstance(node_color, pd.DataFrame) and not np.array_equal(np.sort(node_color.iloc[:, 0].astype(str) + SEP1 + node_color.iloc[:, 1].astype(str)), np.sort(expanded_label[~blank_flag])):
        raise ValueError(
                f"🛑 Please provide a comprehensive combination of datasets and cell types in `node_color`")
    color_mapping = _new_relation_to_color(new_relation, node_color, novel_node_color, remain_node_color)
    node_color = np.array([color_mapping[x] for x in np.where(blank_flag, original_label, expanded_label)], dtype = object)
    #annotations
    datasets = new_relation.columns
    if (len(trace_dict) >= 1) and ('orientation' in trace_dict) and (trace_dict['orientation'] == 'v'):
        annotations = [dict(text = datasets[i], y = 1 - i/(len(datasets)-1), x = 0, xanchor = 'right', xref = 'paper', yanchor = 'top' if i == 0 else ('bottom' if i == len(datasets)-1 else 'middle'), yref = 'paper', showarrow = False) for i in range(len(datasets))]
    else:
        annotations = [dict(text = datasets[i], x = i/(len(datasets)-1), y = 0, yanchor = 'top', yref = 'paper', xanchor = 'left' if i == 0 else ('right' if i == len(datasets)-1 else 'center'), xref = 'paper', showarrow = False) for i in range(len(datasets))]
    #update trace and layout
    fig = go.Figure(trace)
    fig.update_traces(node_label = np.where(blank_flag, '', original_label) if not expand_label else expanded_label, node_color = node_color, link_color = link_color, **trace_dict)
    fig.update_layout(title_text = title, width = width, height = height, annotations = annotations, **layout_dict)
    if show:
        fig.show()
    if save:
        ext = os.path.splitext(save)[1] if isinstance(save, str) else '.html'
        if ext not in ['.html', '.png', '.jpg', '.jpeg', '.webp', '.svg', '.pdf', '.eps']:
            raise ValueError(
                    f"🛑 Please provide valid figure suffix: .html, .png, .jpg, .jpeg, .webp, .svg, .pdf, .eps")
        if ext == '.html':
            fig.write_html(save) if isinstance(save, str) else fig.write_html('CellTypist_sankeyplot.html')
        else:
            fig.write_image(save)

def treeplot(alignment, group_celltype: bool = True, order_dataset: bool = False,
        #link
        link_color: str = '#0000007B', link_width: Optional[float] = None,
        #root and node
        node_shape: Union[list, str] = 'o', node_color: Optional[pd.DataFrame] = None, cmap: Union[matplotlib.colors.Colormap, str] = 'Reds', node_size: Optional[float] = None,
        #label
        show_label: bool = True, label_color: str = '#000000', label_size: Optional[Union[float, str]] = None, label_ha: str = 'center', label_va: str = 'top',
        #figure elements
        title: str = 'CellTypist label harmonization tree',
        #show and/or save figure
        ax: Optional[matplotlib.axes.Axes] = None, figsize: Optional[Union[list, tuple]] = None, show: bool = True, save: Union[str, bool] = False,
        #link setting
        link_dict: dict = {},
        #node and root setting
        node_dict: dict = {},
        #label setting
        label_dict: dict = {},
        #for developer use
        expand_label: bool = False,
        ) -> None:
    """
    Generate a tree showing the CellTypist label harmonization in a qualitative manner.

    Parameters
    ----------
    alignment
        A :class:`~celltypist.contro.align.DistanceAlignment` or :class:`~pandas.DataFrame` object representing the harmonization result.
    group_celltype
        Whether to group cell types (rows) in the harmonization table for plotting.
        (N.B. Do not change the default value of this argument unless you know what you are doing.)
        (Default: `True`)
    order_dataset
        Whether to change the dataset order in each cell type group to manifest as hierarchy (tree).
        (Default: `False`)
    link_color
        Color of links/branches.
        (Default: `'#0000007B'`)
    link_width
        Width of links/branches in points.
        Default to 1.5 in a canonical Matplotlib setting.
    node_shape
        Shape of the node. This can also be a list of symbols for datasets that are aligned.
        (Default: `'o'`)
    node_color
        A :class:`~pandas.DataFrame` with three consecutive columns representing dataset, cell type, and color, respectively.
        Default to a color scheme that allows matched cell types to have the same colors.
        This can also be a data frame with columns of dataset, cell type, and numeric value (for mapping color gradient).
    cmap
        Color map to use. This parameter is only relevant if `node_color` is a value-mapping data frame.
        (Default: `'Reds'`)
    node_size
        Size of nodes (cell types) in points.
        Default to 6.0 in a canonical Matplotlib setting.
    show_label
        Whether to label each node with its cell type name.
        (Default: `True`)
    label_color
        Color of cell type labels.
        (Default: `'#000000'`)
    label_size
        Size of cell type labels.
        Default to 10.0 in a canonical Matplotlib setting.
    label_ha
        Horizontal alignment of cell type labels relative to the nodes.
        (Default: `'center'`)
    label_va
        Vertical alignment of cell type labels relative to the nodes.
        (Default: `'top'`)
    title
        Figure title.
        (Default: `'CellTypist label harmonization tree'`)
    ax
        An :class:`~matplotlib.axes.Axes` where the tree will be drawn. Default to draw the tree on a new axes.
    figsize
        Tuple of figure width and height in inches.
        Default to auto-adjusting the figure size based on the numbers of datasets and cell types.
    show
        Whether to show the plot.
        (Default: `True`)
    save
        Whether to save the plot. This can also be a figure filename.
        (Default: `False`)
    link_dict
        A dict passed to :class:`~matplotlib.lines.Line2D` for setting the links/branches.
    node_dict
        A dict passed to :class:`~matplotlib.lines.Line2D` for setting the nodes.
    label_dict
        A dict passed to :class:`~matplotlib.text.Text` for setting cell type labels.
    expand_label
        Ignored. Whether to show the unique expanded labels. Only for developer use.
        (Default: `False`)

    Returns
    ----------
    None
    """
    if isinstance(alignment, pd.DataFrame):
        relation = alignment
    elif hasattr(alignment, 'relation'):
        relation = alignment.relation
    else:
        raise TypeError(
                f"🛑 Please provide correct input - either a DistanceAlignment or a data frame")
    new_relation = _identify_relation_groups(relation, group_prefix = 'Group', order_row = group_celltype, order_column = order_dataset)[1]
    n_col = new_relation.shape[1]
    n_row = new_relation.shape[0]
    #node coordinates
    node_coord = {}
    for col_rank in range(1, n_col + 1):
        col_content = new_relation.iloc[:, col_rank - 1]
        for celltype in np.unique(col_content):
            rows = n_row - np.where(col_content == celltype)[0]
            node_coord[celltype] = [col_rank, rows.mean()]
    #link pairs
    link_pairs = np.row_stack([new_relation.iloc[:, [i, i+1]].drop_duplicates().values for i in range(n_col - 1)])
    #ax
    if ax is None:
        figsize = [3.5*n_col, n_row/3.5] if figsize is None else figsize
        _, ax = plt.subplots(figsize = figsize)
    #links
    for start, end in link_pairs:
        xs = [node_coord[start][0], node_coord[end][0]]
        ys = [node_coord[start][1], node_coord[end][1]]
        ax.plot(xs, ys, ls = '-', marker = 'None', color = link_color, lw = link_width, **link_dict)
    if order_dataset:
        for first_cell_type in np.unique(new_relation.iloc[:, 0]):
            ax.plot([0.5, 1], [(n_row + 1)/2, node_coord[first_cell_type][1]], ls = '-', marker = 'None', color = link_color, lw = link_width, **link_dict)
    #nodes and labels
    expanded_labels = np.array(list(node_coord.keys()))
    original_labels = pd.Series(expanded_labels).str.split(SEP1, expand = True)[1].values
    blank_flag = np.isin(original_labels, [NOVEL, REMAIN])
    if isinstance(node_color, pd.DataFrame) and not np.array_equal(np.sort(node_color.iloc[:, 0].astype(str) + SEP1 + node_color.iloc[:, 1].astype(str)), np.sort(expanded_labels[~blank_flag])):
        raise ValueError(
                f"🛑 Please provide a comprehensive combination of datasets and cell types in `node_color`")
    color_mapping = _new_relation_to_color(new_relation, node_color, None, None, cmap)
    if not isinstance(node_shape, str) and len(node_shape) < n_col:
        raise ValueError(
                f"🛑 Please provide `node_shape` of length {n_col}")
    node_shapes = dict(zip(relation.columns[::2], [node_shape]*n_col)) if isinstance(node_shape, str) else dict(zip(relation.columns[::2], node_shape))
    for node, coord in node_coord.items():
        dst = node.split(SEP1)[0]
        original_label = node.split(SEP1)[1]
        if original_label in [NOVEL, REMAIN]:
            continue
        ax.plot(coord[0], coord[1], marker = node_shapes[dst], ms = node_size, color = color_mapping[node], ls = 'None', **node_dict)
        if show_label:
            ax.text(coord[0], coord[1], node if expand_label else original_label, color = label_color, size = label_size, ha = label_ha, va = label_va, **label_dict)
    if order_dataset:
        ax.plot(0.5, (n_row + 1)/2, marker = 'o', ms = node_size, color = '#000000', ls = 'None', **node_dict)
    #others
    ax.set(xlim = [0, n_col+1], ylim = [0, n_row+1], title = title)
    ax.set_axis_off()
    if not order_dataset:
        for col_rank in range(1, n_col + 1):
            ax.text(col_rank, n_row+0.5, new_relation.columns[col_rank-1], color = label_color, size = label_size, ha = 'center', va = 'center', weight = 'bold')
    #show and save
    if save:
        plt.savefig(save) if isinstance(save, str) else plt.savefig('CellTypist_treeplot.pdf')
    if show:
        plt.show()
    if save:
        plt.close()

def heatmap(alignment, plot_type: str = 'similarity',
        #colors
        dataset_color: Optional[pd.DataFrame] = None, celltype_color: Optional[pd.DataFrame] = None,
        #cell
        vmin: Optional[float] = None, vmax: Optional[float] = None, cmap: str = 'RdBu_r',
        #labels
        dataset_celltype_sep: str = ': ',
        #dendrogram
        cluster: bool = True, show_row_dendrogram: bool = True, show_col_dendrogram: bool = True,
        #figure
        figsize: Union[list, tuple] = (10, 10), ax: Optional[matplotlib.axes.Axes] = None, show: bool = True, save: Union[str, bool] = False,
        #others
        **kwargs) -> None:
    """
    Generate a heatmap showing the cell type relationships within and across datasets (i.e., meta-analysis).

    Parameters
    ----------
    alignment
        A :class:`~celltypist.contro.align.DistanceAlignment` object containing the meta-analysis result.
    plot_type
        The type of heatmap to show, being either cross-dataset cell type transcriptome similarities (`'similarity'`) or membership (`'membership'`).
        (Default: `'similarity'`)
    dataset_color
        A :class:`~pandas.DataFrame` with two consecutive columns representing dataset and color, respectively.
        Default to a CellTypist color cycle.
    celltype_color
        A :class:`~pandas.DataFrame` with three consecutive columns representing dataset, cell type, and color, respectively.
        Default to a color scheme that allows matched cell types to have the same colors.
    vmin
        Minimal value to anchor the colormap.
        Default to 0 unless `plot_type = 'similarity'` and normalization is not performed during cell type harmonization.
    vmax
        Maximal value to anchor the colormap.
        Default to 1 unless `plot_type = 'similarity'` and normalization is not performed during cell type harmonization.
    cmap
        Mapping from data values to color space.
        (Default: `'RdBu_r'`)
    dataset_celltype_sep
        Separator to connect names of data sets and cell types for displaying.
    cluster
        Whether to cluster the rows and columns of the heatmap.
        (Default: `True`)
    show_{row,col}_dendrogram
        Whether to show the row/column dendrogram.
        (Default: `True`)
    figsize
        Tuple of figure width and height in inches.
        Default to 10 inches in both dimensions.
    ax
        An :class:`~matplotlib.axes.Axes` where the heatmap will be drawn. Default to draw the tree on a new axes.
    show
        Whether to show the plot.
        (Default: `True`)
    save
        Whether to save the plot. This can also be a figure filename.
        (Default: `False`)
    others
        All other parameters are the same as :func:`seaborn.clustermap` with selected tags and customized defaults.

    Returns
    ----------
    None
    """
    pass
