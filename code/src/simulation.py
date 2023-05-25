import os

import pandas as pd
from mewpy.simulation import get_simulator
from cell_lines import CellLineIds

class Simulation:
    ''' class with methods to do simulations on GEMs and ecGEMs'''
    @staticmethod
    def fba_ecGEMs(md, exch_const, internal_const, obj_id, maxm, feas, unfeas, all_flx_d, md_nm, flx_pror_rules, **kwargs):
        '''
        - add flux proportion constraints for DNA methylation reactions (described in literature) to human1.12 ecGEM models
        - add internal and media constraints
        - do FBA that minimizes the usage reaction of 'prot_pool'.
        :param md: cobrapy model
        :param exch_const: medium constraints of type {'reaction_id}: (lb, ub)}
        :param internal_const: constraints to internal reactions of type {'reaction_id}: (lb, ub)}
        :param obj_id: dict {'reaction_id': reaction_weight} with what to use as objective
        :param maxm: bool where 'True' is to maximize the objective, otherwise is minimize
        :param feas: number of feasible models from previous iterations
        :param unfeas: number of unfeasible models from previous iterations
        :param all_flx_d: dict with {'model/cell_line_id': {reaction_id: flux_value}} from previous iterations
        :param md_nm: name of model/cell id
        :param flx_pror_rules: bool indicating to apply flux rules or not
        # :param kwargs - flx_rules: list with flux constraint rules, like for e.g. {'name': 'constA', 'flux_ratios': {'prodDNA5hmC': 1, 'prodDNA5fC': -5.6}}
        :param kwargs - enz_pfba: True/False, True when wanting to do an enzymatic pFBA (maximize objective, set the max objective value and minimize total protein)
        :return - all_flx_d: updated dict with {'model/cell_line_id': {reaction_id: flux_value}}
                - feas: updated number of feasible models
                - None: when infeasible
        '''
        if flx_pror_rules:
            mdreac = [r.id for r in md.reactions]
            constList = list()
            for id in ['arm_prodDNA5hmC', 'prodDNA5hmCNo1']:
                if id in mdreac:
                    flx_prodDNA5hmC = md.reactions.get_by_id(id).flux_expression
                    break
            for id in ['arm_prodDNA5CaC', 'prodDNA5CaCNo1']:
                if (id in mdreac):
                    flx_prodDNA5CaC = md.reactions.get_by_id(id).flux_expression
                    constB = md.problem.Constraint(flx_prodDNA5hmC - (10.2 * flx_prodDNA5CaC), lb=0, ub=0)
                    constList.append(constB)
                    break
            for id in ['arm_prodAPsite4', 'prodAPsite4No1']:
                if id in mdreac:
                    flx_prodAPsite4 = md.reactions.get_by_id(id).flux_expression
                    constC = md.problem.Constraint(flx_prodAPsite4 - ((2 / 3) * flx_prodDNA5CaC), lb=0, ub=0)
                    constList.append(constC)
                    break
            md.add_cons_vars(constList)
        sim = get_simulator(model=md, envcond=exch_const)
        sim.objective = obj_id
        flx = {}
        res = sim.simulate(method='FBA', maximize=maxm, constraints=internal_const)
        if res.status.name == 'OPTIMAL':
            if 'enz_pfba' in kwargs:
                bm = res.fluxes['adaptbiomass']
                internal_const.update({'adaptbiomass': (float(bm), float(bm))})
                sim.objective = {'prot_pool_exchange': 1.0}
                flx = {}
                res = sim.simulate(method='FBA', maximize=False, constraints=internal_const)
            flx = res.fluxes
            flx = {k: 0 if v < 1e-9 else v for k, v in flx.items()}
            feas += 1
            all_flx_d[md_nm] = flx
            return feas, all_flx_d
        else:
            if flx_pror_rules:
                md.remove_cons_vars(constList)
                sim = get_simulator(model=md, envcond=exch_const)
                sim.objective = obj_id
                flx = {}
                res = sim.simulate(method='FBA', maximize=maxm, constraints=internal_const)
                if res.status.name == 'OPTIMAL':
                    flx = res.fluxes
                    flx = {k: 0 if v < 1e-9 else v for k, v in flx.items()}
                    feas += 1
                    all_flx_d[md_nm] = flx
                    return feas, all_flx_d
                else:
                    print(md_nm, 'infeasible')
                    unfeas += 1
                    return unfeas
            else:
                print(md_nm, 'infeasible')
                unfeas += 1
                return unfeas

    # @staticmethod
    # def minSumFluxes_add_rules_GEMs(model, envcond, constraints):
    #     '''
    #     - uses mewpy to do a pFBA without a second objective (i.e. to minimize the sum of fluxes of all reactions withOUT maximizing biomass)
    #     :param model: cobrapy loaded model
    #     :param envcond: dictionary with exchange reactions ids and corresponding bounds
    #     :param constraints: dictionary with internal reactions ids and corresponding bounds
    #     :return: solution instance containing the resulting pfba fluxes and status
    #     '''
    #     sim = get_simulator(model=model, envcond=envcond)
    #     solver = solver_instance(sim)
    #     for r_id in sim.reactions:
    #         lb, _ = sim.get_reaction_bounds(r_id)
    #         if lb < 0:
    #             pos, neg = r_id + '+', r_id + '-'
    #             solver.add_variable(pos, 0, inf, update=False)
    #             solver.add_variable(neg, 0, inf, update=False)
    #     solver.update()
    #     for r_id in sim.reactions:
    #         lb, _ = sim.get_reaction_bounds(r_id)
    #         if lb < 0:
    #             pos, neg = r_id + '+', r_id + '-'
    #             solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0, update=False)
    #             solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0, update=False)
    #     solver.update()
    #     # add flux constraints' rules:
    #     solver.add_constraint('constA', {'prodDNA5hmC': 1, 'prodDNA5fC': -5.6}, '=', 0, update=False)
    #     if 'prodDNA5CaC' in sim.reactions:
    #         solver.add_constraint('constB', {'prodDNA5hmc': 1, 'prodDNA5CaC': -10.2}, '=', 0, update=False)
    #     if 'prodAPsite4' in sim.reactions:
    #         solver.add_constraint('constC', {'prodAPsite4': 1, 'prodDNA5CaC': -(2/3)}, '=', 0, update=False)
    #     solver.update()
    #     objective = dict()
    #     for r_id in sim.reactions:
    #         lb, _ = sim.get_reaction_bounds(r_id)
    #         if lb < 0:
    #             pos, neg = r_id + '+', r_id + '-'
    #             objective[pos] = 1
    #             objective[neg] = 1
    #         else:
    #             objective[r_id] = 1
    #     solution = solver.solve(objective, minimize=True, constraints=constraints)
    #     return solution

    @staticmethod
    def get_env_const_ecGEMs(m, m_nm, constr_ex, med_ex_ids, **kwargs):
        '''
        - get environmental constraints to apply to ecModel
        :param m: model
        :param m_nm: model/cell id
        :param constr_ex: exchange reactions to which apply flux bound constraints OR False when no env constraints besides opening medium are applied
        :param med_ex_ids: ids of medium exchange reactions
        :param kwargs - envcond: path to .xls file with environmental constraints
                      - sheet_name: name of sheet to use from .xls file with environmental constraints
                      - smp_info: path to file with cell line info
        :return - exch_const: dict {'reaction_id': (lb, ub)} with constraints for exchange reactions
                - constr_ex_name: string indicating type of constraint applied ('medium only' or group of reactions where numerical constraints were applied)
        '''
        if constr_ex:
            envcond_df = pd.read_excel(kwargs['envcond'], sheet_name=kwargs['sheet_name'], index_col=2)
            envcond_df = envcond_df.dropna()
            CellLineIds.convert_cell_line_names(file=envcond_df, smp_info=kwargs['smp_info'])
            envcond_df = envcond_df.iloc[:, 1:]
            envcond_df = envcond_df.loc[constr_ex, :]
            constr_ex_name = '-'.join(constr_ex)
            exch_const = dict()
            for con_ex in constr_ex:
                lb_s = envcond_df.loc[con_ex, m_nm]
                ub_s = envcond_df.loc[con_ex, m_nm + '.1']
                if lb_s > ub_s:  # sometimes tables report bounds not in a specific order, so column ending in '.1' might actually report sometimes the lb
                    lb = ub_s
                    ub = lb_s
                else:
                    lb = lb_s
                    ub = ub_s
                if lb < 0 and ub <= 0:
                    exch_const[con_ex + '_REV'] = (-1 * ub, -1 * lb)
                    exch_const[con_ex] = (0.0, 0.0)
                elif lb >= 0 and ub > 0:
                    exch_const[con_ex] = (lb, ub)
                    exch_const[con_ex + '_REV'] = (0.0, 0.0)
                else:
                    exch_const[con_ex + '_REV'] = (0, -1 * lb)
                    exch_const[con_ex] = (0, ub)
            for ex in m.exchanges:
                if ex.id not in exch_const:
                    if ex.id.endswith('REV') and (ex.id.split('_')[0] in med_ex_ids):
                        exch_const[ex.id] = (0.0, float('inf'))
                    elif ex.id.endswith('REV') and (ex.id.split('_')[0] not in med_ex_ids):
                        exch_const[ex.id] = (0.0, 0.0)
                    else:
                        exch_const[ex.id] = (0.0, float('inf'))
            exch_const = {ex: v for ex, v in exch_const.items() if ex in {r.id for r in m.reactions}}
        elif not constr_ex:
            # open entry of medium metabolites and close the entry of the remaining exchange reactions - by default split reversible reactions are (0, inf):
            exch_const = {ex.id: (0.0, 0.0) for ex in m.exchanges if (ex.id.endswith('_REV') and (ex.id.split('_')[0] not in med_ex_ids))}
            constr_ex_name = 'medonly'  # means medium only
        return exch_const, constr_ex_name

    @staticmethod
    def run_FBA(simul, gr_const, maxm, internal_const, feas, unfeas, all_flx_d, md_nm):
        '''
        - runs FBA with mewpy
        :param simul: mewpy simul object
        :param gr_const: whether to constraint cell growth with experimental growth rates
        :param maxm: whether to maximize or minimize the objective
        :param internal_const: dict {reaction_id: (lb, ub)} with internal constraints to apply
        :param feas: number of feasible models from previous iterations
        :param unfeas: number of unfeasible models from previous iterations
        :param all_flx_d: dict with {'model/cell_line_id': {reaction_id: flux_value}} from previous iterations
        :param md_nm: name of model/cell id
        :return - all_flx_d: updated dict with {'model/cell_line_id': {reaction_id: flux_value}}
                - feas: updated number of feasible models
                - None: when infeasible
        '''

        flx = {}  # to guarantee that if 'flx' variable is not obtained in one iteration (infeasible) the value of previous iteration is used
        if gr_const:
            sol = simul.simulate(method='FBA', maximize=maxm, constraints=internal_const)
        else:
            sol = simul.simulate(method='FBA', maximize=maxm)
        flx = sol.fluxes
        flx = {k: 0 if v < 1e-9 else v for k, v in flx.items()}
        if sol.status.name == 'OPTIMAL':
            feas += 1
            all_flx_d[md_nm] = flx
            return feas, all_flx_d
        else:
            print(md_nm, 'infeasible')
            unfeas +=1

    @staticmethod
    def get_exc_grw_flx(m, m_nm, exp_df, all_flx_d, exc_flx, mass_final_flx, md_type):
        '''
        - update dictionaries with simulated exchange reactions' flux values or growth rates
        :param m: model
        :param m_nm: model/cell line id
        :param exp_df: dataframe with experimental values for exchange fluxes. rows are exchange reactions ids
        :param all_flx_d: dict {model_id: {reaction_id: flux value}} from previous iterations
        :param exc_flx: dict {model_id: {exchange_id: flux value}} from previous iterations
        :param mass_final_flx: dict {model_id: {biomass_id: flux value}} from previous iterations
        :param md_type: either 'ecGEM' or 'GEM'
        :return - exc_flx: updated dict (see above)
                - mass_final_flx: updated dict (see above)
        '''
        inD = dict()
        massD = dict()
        mdreac = [r.id for r in m.reactions]
        if md_type == 'ecGEM':
            for rid in exp_df.index:
                if ((rid in mdreac) and (rid + '_REV' in mdreac)):
                    inD[rid] = all_flx_d[m_nm][rid] - all_flx_d[m_nm][rid + '_REV']
                elif rid in mdreac:
                    inD[rid] = all_flx_d[m_nm][rid]
                elif rid+'_REV' in mdreac:
                    inD[rid] = all_flx_d[m_nm][rid + '_REV']
                else:
                    inD[rid] = 0.0
        elif md_type == 'GEM':
            for rid in exp_df.index:
                if rid in mdreac:
                    inD[rid] = all_flx_d[m_nm][rid]
                else:
                    inD[rid] = 0.0
        exc_flx[m_nm] = inD
        massD['adaptbiomass'] = all_flx_d[m_nm]['adaptbiomass']
        mass_final_flx[m_nm] = massD
        return exc_flx, mass_final_flx

    @staticmethod
    def save_all_flx_mds(all_flx_d, obj_id, flx_fld, algo, constr_ex_name, with_tsk, constr_ex, gr_const, cl_spec, demeth_tsk=False):
        '''
        - save dataframe with all fluxes of all models
        :param all_flx_d: dict {model_id: {reaction_id: flux value}} for all fluxes of all models
        :param obj_id: id of reaction set as metabolic objective
        :param flx_fld: folder where to save simulation results
        :param algo: either 'fastcore' or 'init'. algorithm with which models were reconstructed
        :param constr_ex_name: string indicating type of constraint applied ('medium only' or group of reactions where numerical constraints were applied)
        :param with_tsk: bool, whether highest reaction score was given to necessary reaction of tissue specific tasks or not
        :param constr_ex: list of exchange reactions to which apply flux bound constraints. If False, no constraints are applied.
        :param gr_const: whether to constraint growth with experimental growth rates
        :param cl_spec: bool, whether to use cell line specific reaction of 'prodDNAtot' (True) or use the generic instead (False).
        :param demeth_tsk: bool, whether necessary reactions of demethylation tasks were included
                       (included only if those are done in specific cell line)
                       when the reactions needed for other cell specific tasks are NOT included.
                       default is False.
        :return: path to folder where to save results
        '''
        final = pd.DataFrame.from_dict(all_flx_d)
        final.fillna(0, inplace=True)  # NAs correspond to reactions that exist in some models, but not in others
        if with_tsk: fld_fn = os.path.join(flx_fld, algo, 'including_tsks')
        else:
            if demeth_tsk: fld_fn = os.path.join(flx_fld, algo, 'notsk_wdemethtsk')
            else: fld_fn = os.path.join(flx_fld, algo, 'no_tsks')
        if constr_ex: fld_fn = os.path.join(fld_fn, 'w_flux_constr')
        else: fld_fn = os.path.join(fld_fn, 'no_flux_constr')
        if gr_const: fld_fn = os.path.join(fld_fn,'biomass_constr')
        else: fld_fn = os.path.join(fld_fn, 'no_biomass_constr')
        if cl_spec: fld_fn = os.path.join(fld_fn, 'cl_spec_DNAtot')
        else: fld_fn = os.path.join(fld_fn, 'generic_DNAtot')
        if obj_id:
            obj_id_name = '_'.join({k + '{:.0e}'.format(v) for k, v in obj_id.items()})
            flx_path = os.path.join(fld_fn, f'{obj_id_name}_{constr_ex_name}_FBA_fluxes.tsv')
        else:
            flx_path = os.path.join(fld_fn, f'minallflux_{constr_ex_name}_FBA_fluxes.tsv')
        # fld_fn = '/'.join(flx_path.split('/')[0:-1])
        if not os.path.exists(fld_fn):
           os.makedirs(fld_fn)
        final.to_csv(flx_path, sep='\t')
        return fld_fn

    @staticmethod
    def get_env_const_GEMs(m, m_nm, constr_ex, med_ex_ids, **kwargs):
        '''
        - get environmental constraints to apply to GEM
        :param m: model
        :param m_nm: model/cell id
        :param constr_ex: exchange reactions to which apply flux bound constraints OR False when no env constraints besides opening medium are applied
        :param med_ex_ids: ids of medium exchange reactions
        :param kwargs - envcond: path to .xls file with environmental constraints
                      - sheet_name: name of sheet to use from .xls file with environmental constraints
                      - smp_info: path to file with cell line info
        :return - exch_const: dict with environmental constraints to use, {reaction_id: (lb, ub)}
                - constr_ex_name: string indicating type of constraint applied ('medium only' or group of reactions where numerical constraints were applied)
        '''
        if constr_ex:
            envcond_df = pd.read_excel(kwargs['envcond'], sheet_name=kwargs['sheet_name'], index_col=2)
            envcond_df = envcond_df.dropna()
            CellLineIds.convert_cell_line_names(file=envcond_df, smp_info=kwargs['smp_info'])
            envcond_df = envcond_df.iloc[:, 1:]
            envcond_df = envcond_df.loc[constr_ex, :]
            constr_ex_name = '-'.join(constr_ex)
            exch_const = dict()
            for con_ex in constr_ex:
                lb_s = envcond_df.loc[con_ex, m_nm]
                ub_s = envcond_df.loc[con_ex, m_nm + '.1']
                if lb_s > ub_s:  # sometimes tables report bounds not in a specific order, so column ending in '.1' might actually report sometimes the lb
                    lb = ub_s
                    ub = lb_s
                else:
                    lb = lb_s
                    ub = ub_s
                exch_const[con_ex] = (lb, ub)
            for ex in m.exchanges:
                if ex.id not in exch_const:
                    if ex.id in med_ex_ids:
                        exch_const[ex.id] = (-1000.0, 1000.0)
                    else:
                        exch_const[ex.id] = (0.0, 1000.0)
            exch_const = {ex: v for ex, v in exch_const.items() if ex in {r.id for r in m.reactions}}
        elif not constr_ex:
            # open entry of medium metabolites and close the entry of the remaining exchange reactions - by default split reversible reactions are (0, inf):
            exch_const = {ex.id: (-1000.0, 1000.0) if ex.id in med_ex_ids else (0.0, 1000.0) for ex in m.exchanges}
            constr_ex_name = 'medonly'  # means medium only
        return exch_const, constr_ex_name





