import copy
import sys, os

home_path = '/home/tbarata'
sys.path.extend([os.path.join(home_path, pkg, 'src') for pkg in ('cobamp', 'troppo')])

from troppo.tasks.core import Task
from troppo.tasks.task_io import JSONTaskIO
import re
import pandas as pd
import numpy as np
from cobra import Reaction
from mewpy.solvers import solver_instance
from mewpy.simulation import get_simulator
from math import inf
from multiprocessing import cpu_count


class MetTask():
    '''
    The MetTask class aims to do operations on tasks
    '''

    def __init__(self, task_file_path, model):
        '''
        :param task_file_path: path to file with original tasks
        :param model: generic model
        '''
        self.model = model
        # model metabolite_name - metabolite_id dict:
        self.mtb_dct = {mt.name: mt.id[:-1] for mt in model.metabolites}
        self.mtb_dct['ALLMETSIN'] = 'ALLMETSIN'
        self.mtb_list_comp = [mt.id for mt in model.metabolites]
        self.tasks = pd.read_csv(task_file_path, sep='\t')
        self.total_processes = int(cpu_count() - (cpu_count() / 4))

    def savetask(self, *args):
        ''' saves the task in troppo format'''
        task = Task(
            should_fail=args[0],
            inflow_dict=args[1],
            outflow_dict=args[2],
            reaction_dict=args[3],
            name=str(args[4]),
            flux_constraints={},
            mandatory_activity=[],
            annotations={'subsystem': args[5],
                         'system': args[6],
                         'id': str(args[4]),
                         'description': args[7]}
        )
        print(task)
        TasksLst = args[8]
        TasksLst.append(task)
        return TasksLst

    def flux_def(self, ins_out, row, d):
        ''' builds inflow_dict or outflow_dict with the format d[metabolite_id] =[lb, ub]'''
        if pd.isna(row[ins_out + ' LB']):
            v_LB = 0.0
        else:
            v_LB = row[ins_out + ' LB']
        if pd.isna(row[ins_out + ' UB']):
            v_UB = 1000.0
        else:
            v_UB = row[ins_out + ' UB']
        splited = row[ins_out].split(';')
        for el in splited:
            mtb_id = self.mtb_dct[el[:-3]] + el[-2]
            d[mtb_id] = [v_LB, v_UB]
        return d

    def EQUreac(self, row, EQUDict, ReactionDic, id):
        '''
        builds reaction_dict for a EQU reactions of a task with format:
        {task1_0: {{mtbA: coef, mtbB: coef, ...}, (lb, ub)},
        task1_1: {{mtbF: coef, mtbH: coef, ...}, (lb, ub)},
        ...}
        '''
        if re.search(r'\s=>\s|\s<=>\s', row['EQU']):
            splited = re.split(r'\s=>\s|\s<=>\s', row['EQU'])
            substrates = re.split(r'\s\+\s', splited[0])
            prod = re.split(r'\s\+\s', splited[1])
            for el in substrates:
                match = re.search(r'^\d+', el)
                if match:
                    mtbid = self.mtb_dct[el[match.end():-3]] + el[-2]
                    coef = match.group()
                    EQUDict[mtbid] = -1.0 * coef
                else:
                    mtbid = self.mtb_dct[el[:-3]] + el[-2]
                    EQUDict[mtbid] = -1.0
            for el in prod:
                match = re.search(r'^\d+', el)
                if match:
                    mtbid = self.mtb_dct[el[match.end():-3]] + el[-2]
                    coef = match.group()
                    EQUDict[mtbid] = coef
                else:
                    mtbid = self.mtb_dct[el[:-3]] + el[-2]
                    EQUDict[mtbid] = 1.0
        if '<=>' in row['EQU']:
            lb_inf = -1000.0
            ub_inf = 1000.0
        elif re.search(r'\s=>\s', row['EQU']):
            lb_inf = 0.0  # when there is nan in the table
            ub_inf = 1000.0
            lb_fin = row['EQU LB']
            ub_fin = row['EQU UB']
        if pd.isna(row['EQU LB']):
            v_LB = lb_inf
        else:
            v_LB = lb_fin
        if pd.isna(row['EQU UB']):
            v_UB = ub_inf
        else:
            v_UB = ub_fin
        tp_bounds = (v_LB, v_UB)
        ReactionDic[str(id) + '_' + str(len(ReactionDic))] = (EQUDict, tp_bounds)
        return ReactionDic

    def parse_task_file(self, sys_info='none'):
        '''
        creates a list of tasks with adequate formats to be processed by troppo package
        :param sys_info: - 'none' there is no (sub)system info is provided for the tasks (default)
                         - 'inside' (sub)system info is provided in the file with tasks - file format should be same as file of essential metabolic tasks
                         - a pandas dataframe with (sub)system info for the tasks. index is task description and column is 'SYSTEM' and 'SUBSYSTEMS'
        '''
        # get rid of tasks with metabolites whose ids were not found in human1 (like NA[c] or NA[r]):
        # ntf_mtb_pos = self.tasks.loc[(self.tasks['IN'].str.contains('NA\[', na=False) | self.tasks['OUT'].str.contains('NA\[', na=False)), 'DESCRIPTION'].index
        ntf_mtb_pos = self.tasks.loc[(self.tasks['IN'].str.startswith('NA[', na=False) | self.tasks['OUT'].str.startswith('NA[', na=False)), 'DESCRIPTION'].index
        # dna = self.tasks.loc[(self.tasks['IN'].str.contains('DNA\[', na=False) | self.tasks['OUT'].str.contains('DNA\[', na=False)), 'DESCRIPTION'].index
        # ntf_mtb_pos = set(ntf_mtb_pos) - dna
        task_r_lst = list()
        for pos in ntf_mtb_pos:
            if pd.isna(self.tasks.loc[pos]['DESCRIPTION']):
                r = pos
                while pd.isna(self.tasks.loc[r]['DESCRIPTION']):
                    r -= 1
                else:
                    task_r_lst.append(self.tasks.loc[r]['DESCRIPTION'])
            else:
                task_r_lst.append(self.tasks.loc[pos]['DESCRIPTION'])
        task_r_set = np.unique(task_r_lst)
        for n in task_r_set:
            print(n)
            iv = list(self.tasks[self.tasks['DESCRIPTION'] == n].index)[0]
            r = iv + 1
            while pd.isna(self.tasks.loc[r]['DESCRIPTION']):
                r += 1
            else:
                fv = r - 1
                self.tasks.drop(self.tasks.index[iv:fv+1], inplace=True)
                self.tasks.reset_index(drop=True, inplace=True)
        # when system or subsystem info is inside the task file - create dict with syst/subsyst initials:
        if type(sys_info) == str and sys_info == 'inside':
            syst_dic = dict()
            for el in self.tasks.iloc[0, 1].split(';'):
                st = el.split(':')
                k = st[0].replace(' ', '', 1)
                v = st[1].replace(' ', '', 1)
                syst_dic[k] = v
            self.tasks = self.tasks.loc[1:, :]
        id = 0
        TasksLst = list()
        for i in range(len(self.tasks)):
            print(i)
            row = self.tasks.iloc[i]
            ## save the previous task if a new.py task will be processed:
            if type(sys_info) == str and sys_info == 'inside':
                lst = [0, 1]
            else:
                lst = [0]
            if i in lst:  # at the beginning there is still no task to save at that point
                pass
            elif pd.notna(row['DESCRIPTION']):
                self.savetask(shd_fail, inflowDic, outflowDic, ReactionDic, id, subs, syst, description, TasksLst)
            ## define system and subsystem:
            # when system or subsystem info is inside the task file:
            if type(sys_info) == str and sys_info == 'inside':
                if row.iloc[0] == '#':
                    subs = row['ID']
                elif pd.isna(row.iloc[0]) and pd.notna(row['ID']):
                    syst = syst_dic[row['ID']]
            # when system or subsystem info is in another file:
            elif type(sys_info) != str and pd.notna(row['DESCRIPTION']):
                syst = sys_info.loc[row['DESCRIPTION'],'SYSTEM']
                subs = sys_info.loc[row['DESCRIPTION'],'SUBSYSTEM']
            # when there is no system or subsystem info to add:
            elif type(sys_info) == str and sys_info == 'none':
                subs = {}
                syst = {}
            ## define should_fail:
            if pd.notna(row['DESCRIPTION']) and (row['SHOULD FAIL'] == True):
                shd_fail = True
            elif pd.notna(row['DESCRIPTION']) and pd.isna(row['SHOULD FAIL']):
                shd_fail = False
            ## define description and id (task number):
            if pd.notna(row['DESCRIPTION']):
                description = row['DESCRIPTION']
                id += 1
            ## define influx metabolites:
            if pd.notna(row['IN']) and pd.notna(row['DESCRIPTION']):
                inflowDic = dict()
                inflowDic = self.flux_def(ins_out='IN', row=row, d=inflowDic)
            elif pd.notna(row['IN']) and pd.isna(row['DESCRIPTION']):
                inflowDic = self.flux_def(ins_out='IN', row=row, d=inflowDic)
            ## define outflux metabolites:
            if pd.notna(row['OUT']) and pd.notna(row['DESCRIPTION']):
                outflowDic = dict()
                outflowDic = self.flux_def(ins_out='OUT', row=row, d=outflowDic)
            elif pd.notna(row['OUT']) and pd.isna(row['DESCRIPTION']):
                outflowDic = self.flux_def(ins_out='OUT', row=row, d=outflowDic)
            ## define other reactions:
            if pd.notna(row['EQU']) and pd.notna(row['DESCRIPTION']):
                ReactionDic = dict()
                EQUDict = dict()
                ReactionDic = self.EQUreac(row, EQUDict, ReactionDic, id)
            elif pd.notna(row['EQU']) and pd.isna(row['DESCRIPTION']):
                EQUDict = dict()
                ReactionDic = self.EQUreac(row, EQUDict, ReactionDic, id)
            elif pd.notna(row['DESCRIPTION']) and pd.isna(row['EQU']):
                ReactionDic = {}
            ## save the last task:
            if i == (len(self.tasks) - 1):
                TasksLst = self.savetask(shd_fail, inflowDic, outflowDic, ReactionDic, id, subs, syst, description, TasksLst)
        ## sometimes all metabolites 'IN' are also on 'OUT' (represented by 'ALLMETSIN'
        ## to deal with that:
        for t in TasksLst:
            if 'ALLMETSINe' in t.outflow_dict.keys():
                for k, v in t.inflow_dict.items():
                    t.outflow_dict[k] = v
                del (t.outflow_dict['ALLMETSINe'])
        return TasksLst

    def apply_task(self,t, model):
        '''
        applies a task to a model with knockout exchange reactions
        :param t: a task
        :param model: model after knockout of exchange reactions
        :return model: model with added task reactions
                rc_add: list of reactions added
        '''
        # ids of added reactions of a task are kept in this list
        rc_add = list(map(lambda x: 'IN_' + x, t.inflow_dict.keys())) + list(map(lambda x: 'OUT_' + x, t.outflow_dict.keys())) + list(map(lambda x: 'EQU_' + x, t.reaction_dict))
        for mt_id, bds in t.inflow_dict.items():
            print(mt_id, (bds))
            rc = Reaction(id='IN_' + mt_id, lower_bound=bds[0], upper_bound=bds[1])
            model.add_reactions([rc])
            model.reactions.get_by_id('IN_' + mt_id).add_metabolites({mt_id: 1.0})
        for mt_id, bds in t.outflow_dict.items():
            model.add_boundary(model.metabolites.get_by_id(mt_id), reaction_id='OUT_' + mt_id, lb=bds[0], ub=bds[1], type='')
        if t.reaction_dict:  # if EQU reaction exists
            for rc_id, v in t.reaction_dict.items():
                rc = Reaction(id='EQU_' + rc_id, lower_bound=v[1][0], upper_bound=v[1][1])
                model.add_reactions([rc])
                model.reactions.get_by_id('EQU_' + rc_id).add_metabolites(v[0])
        return model, rc_add

    @classmethod
    def minSumFluxes(cls, model, envcond, constraints, **kwargs):
        '''
        - uses mewpy to do a pFBA without a second objective (i.e. to minimize the sum of fluxes of all reactions withOUT maximizing biomass)
        - may as well set additional provided flux constraint rules like e.g. fluxA = (2/3)fluxB
        :param model: cobrapy loaded model
        :param envcond: dictionary with exchange reactions ids and corresponding bounds
        :param constraints: dictionary with internal reactions ids and corresponding bounds
        :param kwargs - flx_rules: list with flux constraint rules, like for e.g. {'name': 'constA', 'flux_ratios': {'prodDNA5hmC': 1, 'prodDNA5fC': -5.6}}
        :return: solution instance containing the resulting pfba fluxes (.values) and status (status.name)
        '''
        sim = get_simulator(model=model, envcond=envcond)
        solver = solver_instance(sim)
        # add variable for each direction of each reversible reaction:
        for r_id in sim.reactions:
            lb, _ = sim.get_reaction_bounds(r_id)
            if lb < 0:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_variable(pos, 0, inf, update=False)
                solver.add_variable(neg, 0, inf, update=False)
        solver.update()
        # set the bounds of those variables:
        for r_id in sim.reactions:
            lb, _ = sim.get_reaction_bounds(r_id)
            if lb < 0:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0, update=False)
                solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0, update=False)
        solver.update()
        # add flux constraints' rules if provided:
        if 'flx_rules' in kwargs:
            for rl in kwargs['flx_rules']:
                solver.add_constraint(rl['name'], rl['flux_ratios'], rl['sig'], rl['val'], update=False)
            solver.update()
        # set objective:
        objective = dict()
        for r_id in sim.reactions:
            lb, _ = sim.get_reaction_bounds(r_id)
            if lb < 0:
                pos, neg = r_id + '+', r_id + '-'
                objective[pos] = 1
                objective[neg] = 1
            else:
                objective[r_id] = 1
        solution = solver.solve(objective, minimize=True, constraints=constraints)
        return solution

    def find_reactions_used(self, tsk_list, model):
        '''
        - identify tasks should_fail=False that work on the generic model
          and corresponding reactions needed for that task to be executed
        - identify tasks should_fail=True that work on the generic model (and shouldn't), in case they exist,
          and corresponding model reactions that allow the task to work
        :param tsk_list: list of tasks, each in troppo format
        :param model: generic model
        :param to_keep_should_do: successful tasks that are should_fail=False and aforementioned reactions
        :param should_not_do_that_passed: tasks that worked but shouldn't and aforementioned reactions
        '''
        to_keep_should_do = list()
        should_not_do_that_passed = list()
        # we do not copy the model and knockout the exchanges on the copy cause it was giving weird behavior when applying the tasks on the copy (even with deep.copy):
        init_ex_bds = []
        for ex in model.exchanges: init_ex_bds.append((ex.id, ex.bounds))
        for tsk in tsk_list:
            with model as md:
                print(tsk.name)
                # apply task to generic model:
                tsk_mtb = list(tsk.inflow_dict.keys()) + list(tsk.outflow_dict.keys())
                if tsk.reaction_dict:
                    for v in tsk.reaction_dict.values():
                        tsk_mtb += list(v[0].keys())
                if [mtb for mtb in tsk_mtb if mtb not in self.mtb_list_comp]:
                    continue  # sometimes the model does not have the compartimentalized metabolite
                else:
                    md, rc_added = self.apply_task(t=tsk, model=md)
                    # the environmental conditions corresponding to close the
                    # exchange reactions of generic model (all boundary reactions in human1 are exchange reactions):
                    envcond = {ex.id: (0, 0) for ex in md.exchanges if ex.id not in rc_added}
                    # the constraints on fluxes of reactions of the task:
                    const = {rc: md.reactions.get_by_id(rc).bounds for rc in rc_added}
                # do pFBA without objective to
                # - identify reactions that carry flux when performing the task, when we minimize the sum of fluxes of all reactions in the model:
                # - know if the task can be done in the generic model:
                print('try_b')
                pfba = self.minSumFluxes(model=md, envcond=envcond, constraints=const)
                # there is no objective besides minimize the sum of fluxes of all reactions in the model (sometimes people could also want to maximize biomass)
                print(pfba.message)
                # if task is should_fail = False and it is performed successfully
                # - then keep the task
                # - save (in mandatory_activity) the minimum number of reactions needed for the task to be done:
                if not tsk.should_fail and pfba.message == 'optimal':
                    added_at_flx = [rid for rid in pfba.to_dataframe().index if re.search(r'^EQU_|^IN_|^OUT_', rid)]
                    flx = pfba.to_dataframe().drop(added_at_flx) # drop(rc_added) not enough. cause if EQU is reversible produces EQU+ or EQU- in the flux dataframe
                    tsk_c = copy.deepcopy(tsk)
                    tsk_c.mandatory_activity = list({rid[:-1] if re.search(r'\+$|\-$', rid) else rid for rid in list(flx[abs(flx['value']) > 1E-6].index)}) # contains the reactions that had flux when a task that should be done is successfully done
                    to_keep_should_do.append(tsk_c)
                # if task is should_fail = True and it is performed successfully
                # - save the task in a separate list together with the reactions required to have flux
                #   this indicates which reactions have to be removed from generic model for the task to not happen
                elif tsk.should_fail and pfba.status == 'optimal':
                    flx = pfba.to_dataframe().drop(rc_added)
                    tsk_c = copy.deepcopy(tsk)
                    tsk_c.mandatory_activity = list({rid[:-1] if re.search(r'\+$|\-$', rid) else rid for rid in list(flx[abs(flx['value']) > 1E-6].index) })  # contains the reactions that had flux when a reaction that should fail didn't fail
                    should_not_do_that_passed.append(tsk_c)
                print('try_e')
        if should_not_do_that_passed:
            print(str(len(should_not_do_that_passed)) + ' tasks that should have failed were done')
        else:
            print('All tasks that should have failed, failed')
        if to_keep_should_do:
            initial = len([t for t in tsk_list if not t.should_fail])
            final = str(len(to_keep_should_do))
            print(f' {final} out of {initial} tasks that should have been done, were successful')
        # bring generic model exchange bounds back to normal:
        for id, bd in init_ex_bds:
            model.reactions.get_by_id(id).bounds = bd
        return to_keep_should_do, should_not_do_that_passed

    def save_as_json(self, path, task_list):
        ''' saves list of tasks in json format'''
        return JSONTaskIO().write_task(path, task_list)

    @staticmethod
    def replace_tsk(original, list_rep):
        '''
        - replace tasks in a 'original' list with tasks from another list when the tasks of the former have same description as the later (i.e. the same task)
          but the name of task (task number) in the original is kept
        - useful when there're common tasks between two list (with different reactions) and we know the correct reactions are in list_rep
        :param original: list of tasks where we want some to be replaced for the correct representation of the task
        :param list_rep: list of tasks where all tasks have the correct representation
        :return:
        '''
        list_rep_des = {t.annotations['description']: t for t in list_rep}
        for i in range(len(original)):
            print('iteration:', i)
            ts = original[i]
            dsc = ts.annotations['description']
            if dsc in list_rep_des:
                ts_r = copy.deepcopy(list_rep_des[dsc])
                print(ts_r.name)
                ts_r.name = ts.name  # copied essential task will have the name of the corresponding task in parsed_full
                print(ts_r.name)
                original[i] = ts_r
        return original

    def load_tasks(self, path):
        '''
        if task list is saved in .json it loads it
        :param path: path to file with list of tasks to load
        '''
        if os.path.exists(path):
            return JSONTaskIO().read_task(path)

    def is_mandatory_act_enough(self):
        '''
        - check if the reactions identified as necessary for a task are in fact enough for that task to be done
        :param self.tsks_of_generic_md: list of tasks with mandatory activity complete
        :param self.model: generic model
        '''
        for tsk in self.tsks_of_generic_md:
            with self.model as md:
                print(tsk.name)
                md, rc_added = self.apply_task(t=tsk, model=md)
                # the environmental conditions corresponding to close the
                # exchange reactions of generic model (all boundary reactions in human1 are exchange reactions):
                envcond = {ex.id: (0, 0) for ex in md.exchanges if ex.id not in rc_added}
                # the constraints on fluxes of reactions of the task:
                const = {rc: md.reactions.get_by_id(rc).bounds for rc in rc_added}
                ma = tsk.mandatory_activity
                const.update({r.id: (0, 0) for r in md.reactions if (r.id not in rc_added) or (r.id not in ma)})
                pfba = self.minSumFluxes(model=md, envcond=envcond, constraints=const)
                if pfba.message == 'optimal':
                    print(f'task {tsk.name}: all good')
                else:
                    print(f'task {tsk.name}: mandatory activity is not enough for reaction to occur')

    def tasks_to_keep_reconst(self, path_rc_sc, thr, path_out, tsk_path):
        '''
        - determines which tasks should be protected during model reconstruction for each of the cell types/conditions based on the task metabolic score
          outputs file with 1s (task to be protected) and 0s (task to not be protected). rows are tasks, columns are conditions/cell lines.
        :param path_rc_sc: path to file with reaction scores for each cell type/condition (columns) and reactions (rows)
        :param thr: tasks with metabolic score (MT score) above this threshold are to be protected
        :param path_out: path to the output file
        :param tsk_path: path to list with tasks. each task in troppo format and with mandatory activity
        '''
        tsk_list = self.load_tasks(tsk_path)
        # get reaction scores calculated at 'cl_mds_fastcore_no_tsk.py':
        rsc_df = pd.read_csv(path_rc_sc, sep='\t', index_col=0)
        # calculate the task metabolic score:
        mt_score = pd.DataFrame({cond: {t.name: np.array([rsc_df.loc[rid, cond] for rid in t.mandatory_activity]).mean() for t in tsk_list} for cond in rsc_df})
        # which tasks have metabolic score above the threshold?:
        final = (mt_score > thr).astype(int)
        final.to_csv(path_out, sep='\t')
        return final

    def prepare_tasks(self, path_save, sys_path=None):
        '''
        - gives a message in case tasks should_fail=True are not done by the model,
        - identifies and formats tasks should_fail=False that are done by the model,
          saving also the reactions necessary for each of those tasks to work under 'task.mandatory_activity'
        :param sys_path: path to file with system info
        :param path_save: path to file where to save tasks should_fail=False that are successfully done by the model,
                          which includes the task mandatory activity
        '''
        if sys_path == 'inside':
            # parse tasks file:
            parsed = self.parse_task_file(sys_info='inside')
        else:
            # parse tasks file and add system info on tasks:
            sys_info = pd.read_excel(sys_path, sheet_name='Sheet1')
            sys_info_f = sys_info.dropna(subset=['DESCRIPTION'])[['DESCRIPTION', 'SYSTEM', 'SUBSYSTEM']]
            sys_info_f.set_index('DESCRIPTION', inplace=True)
            parsed = self.parse_task_file(sys_info=sys_info_f)

        # - identify tasks should_fail=False that work on the generic model
        # and corresponding reactions needed for that task to be executed
        # - identify tasks should_fail=True that work on the generic model (and shouldn't), in case they exist,
        # and corresponding model reactions that allow the task to work
        shd_do, shd_not = self.find_reactions_used(tsk_list=parsed, model=self.model) # All tasks that should have failed, failed
        # give a consecutive numbering to tasks names. Note: names in .json tasks file is going to be different from names in original data .txt files:
        lsts = [shd_do, shd_not]
        for lst in lsts:
            for i, tsk in enumerate(lst): tsk.name = str(i+1)
        # save tasks should_fail=False that work in the model:
        fld_pth = '/'.join(path_save.split('/')[:-1])
        if not os.path.exists(fld_pth):
            os.makedirs(fld_pth)
        self.save_as_json(path=path_save, task_list=shd_do)
        self.tsks_of_generic_md = shd_do








