import numpy as np

from keras.preprocessing import sequence
from rdkit import Chem

from rdkit.Chem import Descriptors
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem import rdMolDescriptors


import networkx as nx
from rdkit.Chem import rdmolops

from utils.filter import HashimotoFilter
import shutil,os

from utils import sascorer, SDF2xyzV2
from utils.rdock_test_MP import vinadock_score

smiles_max_len = 82 #MW250:60, MW300:70

def expanded_node(model,state,val,loop_num):

    all_nodes=[]
    position=[]
    position.extend(state)

    get_int_old=[]
    for j in range(len(position)):
        get_int_old.append(val.index(position[j]))

    get_int=get_int_old

    x=np.reshape(get_int,(1,len(get_int)))
    x_pad= sequence.pad_sequences(x, maxlen=smiles_max_len, dtype='int32',
        padding='post', truncating='pre', value=0.)

    for i in range(loop_num):
        predictions=model.predict(x_pad)
        #print "shape of RNN",predictions.shape


        preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
        preds = np.log(preds) / 1.0
        preds = np.exp(preds) / np.sum(np.exp(preds))


        next_probas = np.random.multinomial(1, preds, 1)
        next_int=np.argmax(next_probas)
        #get_int.append(next_int)
        all_nodes.append(next_int)

    print(all_nodes)
    all_nodes=list(set(all_nodes))

    print(all_nodes)
    return all_nodes


def node_to_add(all_nodes,val):
    added_nodes=[]
    for i in range(len(all_nodes)):
        added_nodes.append(val[all_nodes[i]])

    print(added_nodes)

    return added_nodes



def chem_kn_simulation(model,state,val,added_nodes):
    all_posible=[]

    end="\n"
    #val2=['C', '(', ')', 'c', '1', '2', 'o', '=', 'O', 'N', '3', 'F', '[C@@H]', 'n', '-', '#', 'S', 'Cl', '[O-]', '[C@H]', '[NH+]', '[C@]', 's', 'Br', '/', '[nH]', '[NH3+]', '4', '[NH2+]', '[C@@]', '[N+]', '[nH+]', '\\', '[S@]', '5', '[N-]', '[n+]', '[S@@]', '[S-]', '6', '7', 'I', '[n-]', 'P', '[OH+]', '[NH-]', '[P@@H]', '[P@@]', '[PH2]', '[P@]', '[P+]', '[S+]', '[o+]', '[CH2-]', '[CH-]', '[SH+]', '[O+]', '[s+]', '[PH+]', '[PH]', '8', '[S@@+]']
    for i in range(len(added_nodes)):

        position=[]
        position.extend(state)
        position.append(added_nodes[i])
        #print state
        #print position
        #print len(val2)
        total_generated=[]
        new_compound=[]
        get_int_old=[]
        for j in range(len(position)):
            get_int_old.append(val.index(position[j]))
        #当前层数
        get_int=get_int_old

        x=np.reshape(get_int,(1,len(get_int)))
        x_pad= sequence.pad_sequences(x, maxlen=smiles_max_len, dtype='int32',
            padding='post', truncating='pre', value=0.)
        while not get_int[-1] == val.index(end):
            predictions=model.predict(x_pad)
            #print "shape of RNN",predictions.shape
            preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
            preds = np.log(preds) / 1.0
            preds = np.exp(preds) / np.sum(np.exp(preds))
            next_probas = np.random.multinomial(1, preds, 1)
            #print predictions[0][len(get_int)-1]
            #print "next probas",next_probas
            #next_int=np.argmax(predictions[0][len(get_int)-1])
            next_int=np.argmax(next_probas)
            a=predictions[0][len(get_int)-1]
            next_int_test=sorted(range(len(a)), key=lambda i: a[i])[-10:]
            get_int.append(next_int)
            x=np.reshape(get_int,(1,len(get_int)))
            x_pad = sequence.pad_sequences(x, maxlen=smiles_max_len, dtype='int32',
                padding='post', truncating='pre', value=0.)
            if len(get_int)>smiles_max_len:
                break
        total_generated.append(get_int)
        all_posible.extend(total_generated)


    return all_posible



def predict_smile(all_posible,val):


    new_compound=[]
    for i in range(len(all_posible)):
        total_generated=all_posible[i]

        generate_smile=[]

        for j in range(len(total_generated)-1):
            generate_smile.append(val[total_generated[j]])
        generate_smile.remove("&")
        new_compound.append(generate_smile)

    return new_compound


def make_input_smile(generate_smile):
    new_compound=[]
    for i in range(len(generate_smile)):
        middle=[]
        for j in range(len(generate_smile[i])):
            middle.append(generate_smile[i][j])
        com=''.join(middle)
        new_compound.append(com)
        new_compound.append(com)
    #print len(new_compound)

    return new_compound



def check_node_type(new_compound, generated_dict, sa_threshold = 10, rule = 0, radical = False, hashimoto_filter=False, dict_id=1, trial = 1):
    node_index=[]
    valid_compound=[]
    score=[]
    f_list = open('list_docking_pose_%s.txt' % trial, 'a')
    for i in range(len(new_compound)):
        #check dictionary
        #print('check dictionary', 'comp:', new_compound[i], 'check:', new_compound[i] in generated_dict)
        #去掉重复分子
        if new_compound[i] in generated_dict:
            node_index.append(i)
            valid_compound.append(new_compound[i])
            score.append(generated_dict[new_compound[i]])
            print('duplication!!')
            continue

        ko = Chem.MolFromSmiles(new_compound[i])
        if ko!=None:
            # check hashimoto_filter
            if hashimoto_filter:
                hashifilter = HashimotoFilter()
                hf,_ = hashifilter.filter([new_compound[i]])
                print('hashimoto filter check is', hf)
                if hf[0] == 0:
                    continue

            #check SA_score   SA 分数用于表示化合物的合成可达性，数值越低表示化合物越容易合成。
            SA_score = -sascorer.calculateScore(ko)

            #if new_compound[i][-1] == '\n':
            #    continue

            #print('SA_score', SA_score)
            if sa_threshold < -SA_score:
                continue

            #check radical检查化合物中是否存在自由基
            if radical:
                #koh = Chem.AddHs(ko)  ## get ValueError: Sanitization error: Explicit valence for atom # 3 C, 6, is
                try:
                    koh = Chem.AddHs(ko)#添加H原子
                except ValueError:
                    continue

                fw = Chem.SDWriter('radical_check.sdf')
                try:
                    fw.write(koh)
                    fw.close()
                except ValueError:
                    continue
                #从 SDF 文件中读取分子的结构信息。这些信息包括原子坐标、总电荷、自旋多重度等。
                Mol_atom, Mol_CartX, Mol_CartY, Mol_CartZ,TotalCharge, SpinMulti = SDF2xyzV2.Read_sdf('radical_check.sdf')
                print('radical check', SpinMulti)
                #用于检查是否存在自由基。自由基通常具有非常特殊的自旋多重度，因此检查其是否等于2
                if SpinMulti == 2: #2:open
                    continue

            #check Rule of Five
            weight = round(rdMolDescriptors._CalcMolWt(ko), 2)#计算分子的分子量
            #print('weight:',weight)
            logp = Descriptors.MolLogP(ko)#计算分子的分配系数（logP 值）分子在水和油相中的分配程度，常用于描述分子的亲水性和脂溶性
            #print('logp:',logp)
            donor = rdMolDescriptors.CalcNumLipinskiHBD(ko)#评估分子的口服生物利用度
            #print('donor:',donor)
            acceptor = rdMolDescriptors.CalcNumLipinskiHBA(ko)#计算分子的 Lipinski 酸数（HBA 数）
            #print('acceptor:',acceptor)
            rotbonds = rdMolDescriptors.CalcNumRotatableBonds(ko)#计算分子的可旋转键数
            #print('rotbonds:',rotbonds)
            if rule == 1:
                if weight > 500 or logp > 5 or donor > 5 or acceptor > 10:
                    continue
            if rule == 2:
                if weight > 300 or logp > 3 or donor > 3 or acceptor > 3 or rotbonds > 3:
                    continue
            #查找给定分子的环基
            cycle_list = nx.cycle_basis(nx.Graph(rdmolops.GetAdjacencyMatrix(ko)))
            if len(cycle_list) == 0:
                cycle_length =0
            else:
                cycle_length = max([ len(j) for j in cycle_list ])
            if cycle_length <= 6:
                cycle_length = 0
            if cycle_length==0:
                m=vinadock_score(new_compound[i])
                if m[0]<10**10 and m[1]<10**10:
                    node_index.append(i)
                    valid_compound.append(new_compound[i])
                    score.append(m[:2])
                    #print(score)
                    #add dictionary
                    generated_dict[new_compound[i]] = m[:2]
                    #print(generated_dict)
                    compound_id = i
                    f_list.write('pose_'+str(dict_id)+'_'+str(compound_id)+'_'+','+new_compound[i]+','+str(m[0])+','+str(m[1])+','+str(SA_score)+','+str(weight)+','+str(logp)+','+str(donor)+','+str(acceptor)+','+str(rotbonds))
                    f_list.write('\n')
    f_list.close()
    return node_index,score,valid_compound, generated_dict