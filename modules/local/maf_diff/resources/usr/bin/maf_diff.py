#!/usr/bin/env python3

import csv
import argparse
from SecretColors import Palette
import plotly.express as px
import plotly.graph_objects as go

palette = Palette("material")


def load_maf(maf_file):
    maf_dict = {}
    csv_reader = csv.DictReader(filter(lambda row: row[0]!='#', maf_file),delimiter='\t')
    for row in csv_reader:
        position_id = (row['Chromosome'],row['Start_Position'],row['End_Position'])
        maf_dict[position_id] = {
            "type": row['Variant_Type'],
            "class": row['Variant_Classification'],
            "hugo": row['Hugo_Symbol'],
            "allele_R": row['Reference_Allele'],
            "allele_T1": row['Tumor_Seq_Allele1'],
            "allele_T2": row['Tumor_Seq_Allele2']

        }
    return maf_dict

def get_random_color():
    palette.color_mode = "rgba"
    return palette.random()

def create_graph(maf_dicts,label, range):
    source_list = []
    label_list = []
    color_list = []
    link_color = []
    link_label = []
    color_dict = {}
    target_list = []
    value_list = []
    node_dict = {}
    link_dict = {}
    new_node_info = {"num": 0 }
    maf_index = 0
    table_cols = ['Transition', 'Variant Type','From','To','Hugo', 'Chromosome', 'Start', 'End']
    table_data = [[],[],[],[],[],[],[],[]]
    tsv_data = "\t".join(table_cols)
    table_data_items = []
    def create_node(node_label,label, maf_id, hugo, var_type):
        if node_label in node_dict:
            return node_dict[node_label]
        id = new_node_info["num"]
        new_node_info["num"] += 1
        if label in color_dict:
            color = color_dict[label]
        else:
            color = get_random_color()
            color_dict[label] = color
        node_dict[node_label] = {"id":id, "label": label, "maf_id": maf_id, "color": color, "hugo": hugo, "var_type": var_type}
        label_list.append(node_label)
        color_list.append('rgba({},{},{},{})'.format(int(color[0]*255),int(color[1]*255),int(color[2]*255),color[3]))
        return node_dict[node_label]

    def create_link(node1, node2, from_maf, to_maf, chr, start, end):
        if node2['label'] =="Not Found" or node1["label"] != node2["label"]:
            link_id = (node1['id'],node2['id'])
            link_label = "{}->{}".format(node1["label"],node2["label"])
            if link_id not in link_dict:
                node_color = node2["color"]
                new_color = 'rgba({},{},{},{})'.format(int(node_color[0]*255),int(node_color[1]*255),int(node_color[2]*255),.5)
                link_dict[link_id] = {'value':0,'color':new_color,'label':link_label}
            link_dict[link_id]['value'] += 1
            table_item_id = "{}_{}_{}_{}_{}_{}".format(node1['maf_id'], node2['maf_id'], link_label,chr,start,end)
            if table_item_id not in table_data_items:
                table_data_items.append(table_item_id)
                table_data[0].append(link_label)
                table_data[1].append(node1['var_type'])
                table_data[2].append(from_maf)
                table_data[3].append(to_maf)
                table_data[4].append(node1['hugo'])
                table_data[5].append(chr)
                table_data[6].append(start)
                table_data[7].append(end)


    def create_tree_data():
        for link_key,link_obj in link_dict.items():
            source = link_key[0]
            target = link_key[1]
            source_list.append(source)
            target_list.append(target)
            value_list.append(link_obj['value'])
            link_color.append(link_obj['color'])
            link_label.append(link_obj['label'])

    def create_sankey():
        source_value = {}
        index = 0
        while index < len(source_list):
            source = source_list[index]
            value = value_list[index]
            if source in source_value:
                source_value[source] += value
            else:
                source_value[source] = value
            index += 1
        source_ranked_list = []
        for single_source in source_value:
            source_ranked_list.append((single_source,source_value[single_source]))
        source_ranked_list.sort(reverse=True,key=lambda x: x[1])
        top_ranked = []
        for single_ranked_source in source_ranked_list[:10]:
            top_ranked.append(single_ranked_source[0])

        new_source_list = []
        new_target_list = []
        new_value_list = []
        new_link_label = []
        new_link_color = []

        index_link = 0
        while index_link < len(source_list):
            source = source_list[index_link]
            if source in top_ranked:
                new_source_list.append(source)
                new_target_list.append(target_list[index_link])
                new_value_list.append(value_list[index_link])
                new_link_label.append(link_label[index_link])
                new_link_color.append(link_color[index_link])
            index_link += 1

        sankey = go.Figure(data=[go.Sankey(
        valueformat = ".0f",
        valuesuffix = "",
        # Define nodes
        node = dict(
            pad = 15,
            thickness = 15,
            line = dict(color = "black", width = 0.5),
            label =  label_list,
            color =  color_list
            ),
        # Add links
        link = dict(
        source =  new_source_list,
        target =  new_target_list,
        value =  new_value_list,
        label =  new_link_label,
        color =  new_link_color
        ))])

        sankey.update_layout(title_text="Varriant Classifications changes across MAFs",
                        font_size=10)
        return sankey

    def create_table(tsv_data):
        if table_data[0]:
            sorted_table = sorted(zip(table_data[0],table_data[1],table_data[2],table_data[3],table_data[4],table_data[5],table_data[6],table_data[7]))
            col1, col2, col3, col4, col5, col6, col7, col8 = zip(*sorted_table)
        else:
            col1 = col2 = col3 = col4 = col5 = col6 = col7 = col8 = []
        index = 0
        while index < len(col1):
            tsv_data = tsv_data + "\n" + "\t".join([col1[index], col2[index], col3[index], col4[index], col5[index], col6[index], col7[index], col8[index]])
            index += 1
        with open("maf_discordant_events_{}.tsv".format(label), 'w') as table_output:
            table_output.write(tsv_data)
        table = go.Figure(data=[go.Table(header=dict(values=table_cols),cells=dict(values=[col1, col2, col3, col4, col5, col6, col7, col8]))])
        table.update_layout(height=1000)
        return table

    def write_html(sankey_figure,table_figure):
        with open("maf_diff_{}.html".format(label), 'a') as html_output:
            html_output.write(sankey_figure.to_html(full_html=False, include_plotlyjs='cdn'))
            html_output.write(table_figure.to_html(full_html=False, include_plotlyjs='cdn'))

    def find_match(single_pos,single_event,next_maf):
        (chromosome, start, end) = single_pos
        if( single_pos in next_maf and
        next_maf[single_pos]['allele_R'] == single_event['allele_R']  and
        next_maf[single_pos]['allele_T1'] == single_event['allele_T1'] and
        next_maf[single_pos]['allele_T2'] == single_event['allele_T2']):
            next_event = next_maf[single_pos]
            next_node_label = next_maf_id +"_"+ next_event["class"]
            next_label = next_event["class"]
            return create_node(next_node_label,next_label, next_maf_id, next_event['hugo'], next_event['type'])
        start_check = int(start) - range
        start_check_end = int(start) + range
        end_check_end = int(end) + range
        while(start_check < start_check_end):
            end_check = int(end) - range
            while(end_check < end_check_end):
                position = (chromosome,str(start_check), str(end_check))
                if( position in next_maf and
                next_maf[position]['allele_R'] == single_event['allele_R'] and
                next_maf[position]['allele_T1'] == single_event['allele_T1'] and
                next_maf[position]['allele_T2'] == single_event['allele_T2']):
                    next_event = next_maf[position]
                    start_diff = start_check - int(start)
                    end_diff = end_check - int(end)
                    shift_label = ""
                    if start_diff:
                        if start_diff > 5:
                            start_diff = "5+"
                        elif start_diff > 0:
                            start_diff = "<5+"
                        elif start_diff < -5:
                            start_diff = "-5+"
                        else:
                            start_diff = "<5-"
                    if end_diff:
                        if end_diff > 5:
                            end_diff = "5+"
                        elif end_diff > 0:
                            end_diff = "<5+"
                        elif end_diff < -5:
                            end_diff = "-5+"
                        else:
                            end_diff = "<5-"
                    if start_diff and end_diff:
                        shift_label = "start_shift_{}_end_shift_{}".format(start_diff,end_diff)
                    elif end_diff:
                        shift_label = "end_shift_{}".format(end_diff)
                    elif start_diff:
                        shift_label = "start_shift_{}".format(start_diff)
                    next_node_label = next_maf_id +"_"+ next_event["class"] + "_" + shift_label
                    next_label = next_event["class"] + "_" + shift_label
                    return create_node(next_node_label,next_label, next_maf_id, next_event['hugo'], next_event['type'])
                end_check += 1
            start_check += 1
        next_node_label = next_maf_id + "_Not Found"
        next_event = None
        next_label = "Not Found"
        return create_node(next_node_label,next_label, next_maf_id, None, None)






    while (maf_index +1) < len(maf_dicts):
        current_maf_id = maf_dicts[maf_index][1]
        next_maf_id = maf_dicts[maf_index+1][1]
        current_maf = maf_dicts[maf_index][0]
        next_maf = maf_dicts[maf_index+1][0]
        for single_pos,single_event in current_maf.items():
            current_node_label = current_maf_id +"_" +single_event["class"]
            current_label = single_event["class"]
            current_node = create_node(current_node_label,current_label, current_maf_id, single_event['hugo'], single_event['type'])
            next_node = find_match(single_pos,single_event,next_maf)
            create_link(current_node,next_node, current_maf_id, next_maf_id, single_pos[0], single_pos[1], single_pos[2])
        maf_index += 1
    create_tree_data()
    sankey = create_sankey()
    table = create_table(tsv_data)
    write_html(sankey, table)


def diff_mafs(maf_list,label_list,range):
    maf_dicts = []
    index = 0
    label = ""
    while index < len(maf_list):
        single_maf = maf_list[index]
        single_label = label_list[index]
        maf_obj = load_maf(single_maf)
        maf_dicts.append((maf_obj, single_label))
        label = label + "_{}".format(single_label)
        index += 1
    create_graph(maf_dicts,label,range)



if __name__ == '__main__':
    parser = argparse.ArgumentParser("maff_diff")
    parser.add_argument('--mafs',type=argparse.FileType('r', encoding='UTF-8'),help="The maf files", required=True,nargs='+')
    parser.add_argument('--labels',type=str, help="labels for the maf files in the same order", required=True,nargs='+')
    parser.add_argument('--range',type=int, default=10, help="Range to check in both -+ direction for a shifted read", required=False)
    args = parser.parse_args()
    if len(args.mafs) < 2 or len(args.labels) < 2:
        print("At least 2 mafs or labels must be specified")
    if len(args.mafs) != len(args.labels):
        print("Labels and maf lists should be the same size")
    diff_mafs(args.mafs, args.labels, args.range)
