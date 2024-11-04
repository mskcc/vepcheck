#!/usr/bin/env python3

import csv
import argparse
from SecretColors import Palette
import plotly.express as px
import plotly.graph_objects as go
from multiprocessing.pool import ThreadPool
from threading import current_thread, Lock

palette = Palette("material")


def load_maf(maf_file):
    maf_dict = {}
    label_dict = {}
    csv_reader = csv.DictReader(filter(lambda row: row[0]!='#', maf_file),delimiter='\t')
    for row in csv_reader:
        position_id = (row['Chromosome'],row['Start_Position'],row['End_Position'], row['Reference_Allele'], row['Tumor_Seq_Allele1'], row['Tumor_Seq_Allele2'])
        row_id = None
        transcript_id = "NA"
        if "Transcript_ID" in row:
            transcript_id = row["Transcript_ID"]
        if "vcf_id" in row and row["vcf_id"] != ".":
            row_id = row["vcf_id"]
            label_dict[row_id] = position_id
        maf_dict[position_id] = {
            "type": row['Variant_Type'],
            "class": row['Variant_Classification'],
            "hugo": row['Hugo_Symbol'],
            "row_id": row_id,
            "transcript_id": transcript_id
        }
    return (maf_dict, label_dict)

def get_random_color():
    palette.color_mode = "rgba"
    return palette.random()

def create_graph(maf_dicts,label, range, threads, max):
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
    table_cols = ['Transition', 'Variant Type','From','To','Hugo', 'Chromosome', 'Start', 'End', 'Start_shift', 'End_shift']
    stats_dict = {'total': 0, 'shifted': 0, 'match': 0, 'shifted_and_unmatched': 0, 'unmatched': 0, 'missing': 0, 'divergent_transcripts': 0}
    table_data = [[],[],[],[],[],[],[],[],[],[]]
    tsv_data = "\t".join(table_cols)
    table_data_items = []
    node_lock = Lock()
    link_lock = Lock()
    stats_lock = Lock()

    def create_node(node_label,label, maf_id, hugo, var_type, transcript_id):
        if node_label in node_dict:
            return node_dict[node_label]
        id = new_node_info["num"]
        new_node_info["num"] += 1
        if label in color_dict:
            color = color_dict[label]
        else:
            color = get_random_color()
            color_dict[label] = color
        with node_lock:
            node_dict[node_label] = {"id":id, "label": label, "maf_id": maf_id, "transcript_id": transcript_id, "color": color, "hugo": hugo, "var_type": var_type}
            label_list.append(node_label)
            color_list.append('rgba({},{},{},{})'.format(int(color[0]*255),int(color[1]*255),int(color[2]*255),color[3]))
        return node_dict[node_label]

    def create_link(node1, node2, from_maf, to_maf, chr, start, end, start_shift,end_shift):
        with stats_lock:
            if node1['transcript_id'] != node2['transcript_id']:
                stats_dict["divergent_transcripts"] += 1
            if node2['label'] =="Not Found":
                stats_dict["missing"] += 1
            elif node1["label"] != node2["label"]:
                if start_shift != 0 or end_shift != 0:
                    stats_dict["shifted_and_unmatched"] += 1
                else:
                    stats_dict["unmatched"] += 1
            elif start_shift != 0 or end_shift != 0:
                stats_dict["shifted"] += 1
            else:
                stats_dict["match"] += 1
            stats_dict['total'] += 1

        if node1['transcript_id'] == node2['transcript_id']:
            if node2['label'] =="Not Found" or node1["label"] != node2["label"] or start_shift != 0 or end_shift != 0:
                link_id = (node1['id'],node2['id'])
                link_label = "{}->{}".format(node1["label"],node2["label"])
                if link_id not in link_dict:
                    node_color = node2["color"]
                    new_color = 'rgba({},{},{},{})'.format(int(node_color[0]*255),int(node_color[1]*255),int(node_color[2]*255),.5)
                    link_dict[link_id] = {'value':0,'color':new_color,'label':link_label}
                link_dict[link_id]['value'] += 1
                table_item_id = "{}_{}_{}_{}_{}_{}".format(node1['maf_id'], node2['maf_id'], link_label,chr,start,end)
                if table_item_id not in table_data_items:
                    with link_lock:
                        table_data_items.append(table_item_id)
                        table_data[0].append(link_label)
                        table_data[1].append(node1['var_type'])
                        table_data[2].append(from_maf)
                        table_data[3].append(to_maf)
                        table_data[4].append(node1['hugo'])
                        table_data[5].append(chr)
                        table_data[6].append(start)
                        table_data[7].append(end)
                        table_data[8].append(str(start_shift))
                        table_data[9].append(str(end_shift))


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
        for single_ranked_source in source_ranked_list[:max]:
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
            sorted_table = sorted(zip(table_data[0],table_data[1],table_data[2],table_data[3],table_data[4],table_data[5],table_data[6],table_data[7], table_data[8], table_data[9]))
            col1, col2, col3, col4, col5, col6, col7, col8, col9, col10 = zip(*sorted_table)
        else:
            col1 = col2 = col3 = col4 = col5 = col6 = col7 = col8 = col9 = col10 = []
        index = 0
        while index < len(col1):
            tsv_data = tsv_data + "\n" + "\t".join([col1[index], col2[index], col3[index], col4[index], col5[index], col6[index], col7[index], col8[index], str(col9[index]), str(col10[index])])
            index += 1
        with open("maf_discordant_events_{}.tsv".format(label), 'w') as table_output:
            table_output.write(tsv_data)
        table = go.Figure(data=[go.Table(header=dict(values=table_cols),cells=dict(values=[col1, col2, col3, col4, col5, col6, col7, col8, col9, col10]))])
        table.update_layout(height=1000)
        return table

    def write_html(sankey_figure,table_figure):
        total = stats_dict['total']
        matched = stats_dict['match']
        shifted_and_unmatched = stats_dict['shifted_and_unmatched']
        shifted = stats_dict['shifted']
        unmatched = stats_dict['unmatched']
        missing = stats_dict['missing']
        divergent_transcripts = stats_dict['divergent_transcripts']
        stats_line = "<p> Stats: Total events: {}, matched {} ({}%), shifted: {} ({}%), unmatched: {} ({}%), shifted_and_unmatched: {} ({}%), divergent_transcripts: {} ({}%), missing: {} ({}%) </p>".format(str(total), str(matched), str(round(matched/total*100,4)), str(shifted), str(round(shifted/total*100,4)), str(unmatched), str(round(unmatched/total*100,4)), str(shifted_and_unmatched), str(round(shifted_and_unmatched/total*100,4)),str(divergent_transcripts), str(round(divergent_transcripts/total*100,4)), str(missing), str(round(missing/total*100,4)))
        with open("maf_diff_{}.html".format(label), 'a') as html_output:
            html_output.write(stats_line)
            if not table_data[0]:
                if divergent_transcripts == 0:
                    html_output.write("<p> A complete match! </p>")
            else:
                html_output.write(sankey_figure.to_html(full_html=False, include_plotlyjs='cdn'))
                html_output.write(table_figure.to_html(full_html=False, include_plotlyjs='cdn'))

    def find_match(single_pos,single_event,next_maf, next_maf_id, next_maf_row_id_dict):
        maf_row_id = single_event["row_id"]
        (chromosome, row_start, row_end, row_ref, row_allele1, row_allele2) = single_pos
        if maf_row_id and maf_row_id in next_maf_row_id_dict:
            next_pos = next_maf_row_id_dict[maf_row_id]
            next_event = next_maf[next_pos]
            next_label = next_event["class"]
            next_node_label = next_maf_id +"_"+ next_event["class"]
            node = create_node(next_node_label,next_label, next_maf_id, next_event['hugo'], next_event['type'], next_event['transcript_id'])
            start_diff = int(next_pos[1]) - int(row_start)
            end_diff = int(next_pos[2]) - int(row_end)
            return node, start_diff, end_diff
        if single_pos in next_maf:
            next_pos = single_pos
            next_event = next_maf[single_pos]
            next_node_label = next_maf_id +"_"+ next_event["class"]
            next_label = next_event["class"]
            node = create_node(next_node_label,next_label, next_maf_id, next_event['hugo'], next_event['type'], next_event['transcript_id'])
            start_diff = int(next_pos[1]) - int(row_start)
            end_diff = int(next_pos[2]) - int(row_end)
            return node, start_diff, end_diff

        start_check = int(row_start) - range
        start_check_end = int(row_start) + range
        end_check_end = int(row_end) + range
        while(start_check < start_check_end):
            end_check = int(row_end) - range
            while(end_check < end_check_end):
                position = (chromosome,str(start_check), str(end_check), row_ref, row_allele1, row_allele2)
                if position in next_maf:
                    next_event = next_maf[position]
                    start_diff = start_check - int(row_start)
                    end_diff = end_check - int(row_end)
                    next_node_label = next_maf_id +"_"+ next_event["class"] + "_shifted"
                    next_label = next_event["class"] + "_shifted"
                    node = create_node(next_node_label,next_label, next_maf_id, next_event['hugo'], next_event['type'], next_event['transcript_id'])
                    return node, start_diff, end_diff
                end_check += 1
            start_check += 1
        next_node_label = next_maf_id + "_Not Found"
        next_event = None
        next_label = "Not Found"
        node =  create_node(next_node_label,next_label, next_maf_id, 0, 0,0)
        return node, 0, 0

    def process_event(event_list):
        current_maf_id = event_list[0]
        next_maf_id = event_list[2]
        next_maf = event_list[3]
        next_maf_row_id_dict = event_list[4]
        events = event_list[5]
        current_chr = "NA"
        current_thr = current_thread()._name
        if events:
            current_chr = events[0][0][0]
            print("[{}] Working on chromosome {} with {} events".format(current_thr, current_chr, len(events)))

        for single_pos, single_event in events:
            current_node_label = current_maf_id +"_" +single_event["class"]
            current_label = single_event["class"]
            current_node = create_node(current_node_label,current_label, current_maf_id, single_event['hugo'], single_event['type'], single_event['transcript_id'])
            next_node, start_diff, end_diff = find_match(single_pos,single_event,next_maf,next_maf_id,next_maf_row_id_dict)
            create_link(current_node,next_node, current_maf_id, next_maf_id, single_pos[0], single_pos[1], single_pos[2], start_diff, end_diff)
        print("[{}] Finished working on chromosome {}".format(current_thr, current_chr))





    while (maf_index +1) < len(maf_dicts):
        current_maf_id = maf_dicts[maf_index][1]
        current_maf_row_id_dict = maf_dicts[maf_index][2]
        next_maf_id = maf_dicts[maf_index+1][1]
        current_maf = maf_dicts[maf_index][0]
        next_maf = maf_dicts[maf_index+1][0]
        next_maf_row_id_dict = maf_dicts[maf_index+1][2]
        chr_dict = {}
        event_list = []
        for single_pos,single_event in current_maf.items():
            chromosome = single_pos[0]
            if chromosome not in chr_dict:
                chr_dict[chromosome] = []
            chr_dict[chromosome].append((single_pos, single_event))
        for single_chr in chr_dict:
            chr_events = chr_dict[single_chr]
            event_list.append((current_maf_id, current_maf_row_id_dict, next_maf_id, next_maf, next_maf_row_id_dict, chr_events))

        with ThreadPool(threads) as pool:
            pool.map(process_event, event_list)
        maf_index += 1
    create_tree_data()
    sankey = create_sankey()
    table = create_table(tsv_data)
    write_html(sankey, table)


def diff_mafs(maf_list,label_list,range,threads, max):
    maf_dicts = []
    index = 0
    label = ""
    while index < len(maf_list):
        single_maf = maf_list[index]
        single_label = label_list[index]
        maf_obj, maf_row_label = load_maf(single_maf)
        maf_dicts.append((maf_obj, single_label, maf_row_label))
        label = label + "_{}".format(single_label)
        index += 1
    create_graph(maf_dicts,label,range,threads, max)



if __name__ == '__main__':
    parser = argparse.ArgumentParser("maff_diff")
    parser.add_argument('--mafs',type=argparse.FileType('r', encoding='UTF-8'),help="The maf files", required=True,nargs='+')
    parser.add_argument('--labels',type=str, help="labels for the maf files in the same order", required=True,nargs='+')
    parser.add_argument('--range',type=int, default=10, help="Range to check in both -+ direction for a shifted read", required=False)
    parser.add_argument('--threads',type=int, default=10, help="Number of threads to use", required=False)
    parser.add_argument('--max',type=int, default=100, help="Show only a set top amount of classification changes", required=False)
    args = parser.parse_args()
    if len(args.mafs) < 2 or len(args.labels) < 2:
        print("At least 2 mafs or labels must be specified")
    if len(args.mafs) != len(args.labels):
        print("Labels and maf lists should be the same size")
    diff_mafs(args.mafs, args.labels, args.range, args.threads, args.max)
