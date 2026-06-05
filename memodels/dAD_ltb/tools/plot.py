import matplotlib.pyplot as plt

import numpy
import plottools as pt
import bluepyopt
import collections
import seaborn.apionly as sns
import pandas as pd

import logging
logger = logging.getLogger(__name__)

def evolution(log, figs={}, color='b', reportname=''):
    """Plot evolution"""

    figname = 'Evolution ' + reportname
    fig = pt.make_figure(figname=figname, figs=figs, fontsizes = (8,8))

    if figs[figname]['axs'] is not None:
        ax = figs[figname]['axs']
    else:
        ax = pt.single_ax(fig=fig)
        ax.yaxis.grid(which='both')
        figs[figname]['axs'] = ax


    gen_numbers = log.select('gen')
    mean = numpy.array(log.select('avg'))
    std = numpy.array(log.select('std'))
    minimum = numpy.array(log.select('min'))

    stdminus = mean - std
    stdplus = mean + std
    ax.plot(
        gen_numbers,
        mean,
        color=color,
        linewidth=2,
        alpha=0.4,
        label='population average')

    ax.fill_between(
        gen_numbers,
        stdminus,
        stdplus,
        color=color,
        alpha=0.1,
        linewidth=2,
        label=r'population standard deviation')

    ax.plot(
        gen_numbers,
        minimum,
        color=color,
        linewidth=2,
        label='population minimum')

    minmin = min(minimum)
    ax.set_yscale("log")
    ax.set_xlim(min(gen_numbers) - 1, max(gen_numbers) + 1)
    ax.set_xlabel('Generation #')
    ax.set_ylabel('Sum of objectives')
    ax.set_ylim([minmin-0.1*minmin, max(stdplus)])

    ax.axhline(minmin, color='k', linewidth=0.5, label='minimum: %.2f' % minmin)

    ax.legend()

    fig.suptitle(reportname, fontsize=10)

    return fig


def responses(responses, color='b', figs={}, cols=1, reportname='', label=''):
    '''creates subplots'''

    figname = 'Responses ' + reportname

    traces = []
    for response_name, response in responses.iteritems():
        if isinstance(response, bluepyopt.ephys.responses.TimeVoltageResponse):
            traces.append(response_name)

    traces = sorted(traces)
    plot_count = len(traces)

    fig = pt.make_figure(figname=figname,
                            orientation='page',
                            figs=figs,
                            fontsizes = (6,8))

    if figs[figname]['axs'] is not None:
        axs = figs[figname]['axs']
        figs[figname]['i'] += 1
        i = figs[figname]['i']
    else:
        axs = pt.grid_axs(d_out=5, rows=plot_count, columns=cols,
                            fig=fig,
                            top_margin=0.05, bottom_margin=0.05,
                            left_margin=0.12, right_margin=0.05,
                            hspace=0.3, wspace=0.3)
        figs[figname]['axs'] = axs
        i = 0
        figs[figname]['i'] = i

    axs[i].set_title(label)

    for j, response_name in enumerate(traces):

        response = responses[response_name]

        n = j*cols+i
        axs[n].plot(
            response['time'],
            response['voltage'],
            color=color,
            linewidth=1,
            alpha=1)

        axs[n].set_ylabel(response_name.replace('.', '\n'))
        axs[n].set_autoscaley_on(True)
        axs[n].set_autoscalex_on(True)

        #axs[n].set_ylim((-85, 50))
        axs[n].set_xlim((response['time'].as_matrix()[0], response['time'].as_matrix()[-1]))

    axs[-1].set_xlabel('t (ms)')

    fig.suptitle(reportname, fontsize=10)

    return fig


def objectives(objectives, features, color='b', figs={}, reportname='', label=''):
    """Plot objectives and features of the cell model"""

    obj_copy = objectives.copy()
    #del obj_copy['global_maximum']

    max_feature_name = obj_copy.keys()[obj_copy.values().index(max(obj_copy.values()))]
    features['global_maximum'] = features[max_feature_name]

    objectives = collections.OrderedDict(sorted(objectives.iteritems()))
    features = collections.OrderedDict(sorted(features.iteritems()))

    figname = 'Objectives ' + reportname
    fig = pt.make_figure(figname=figname,
                        orientation='page',
                        figs=figs,
                        fontsizes = (8,8))
    
    if figs[figname]['axs'] is not None:
        ax = figs[figname]['axs']
        old_obj_keys = figs[figname]['obj_keys']
    else:
        ax = pt.single_ax(d_out=5, fig=fig,
                    top_margin=0.05, bottom_margin=0.05,
                    left_margin=0.4, right_margin=0.1)
        figs[figname]['axs'] = ax
        old_obj_keys = None

    obj_val = objectives.values()
    obj_keys = objectives.keys()

    feat_val = features.values()
    feat_keys = features.keys()

    
    if old_obj_keys is not None:
        add_obj_keys = list(set(old_obj_keys).symmetric_difference(set(obj_keys)))
        new_obj_keys = old_obj_keys + add_obj_keys
        i_arrange = [obj_keys.index(x) for x in new_obj_keys]
        obj_val = numpy.array(obj_val)[i_arrange]
        obj_keys = new_obj_keys

    ytick_pos = [x + 0.5 for x in range(len(obj_keys))]

    sv = sum(objectives.values())
    logger.info('Sum of objectives %0.4f (# std)', sv)

    ax.barh(ytick_pos,
              obj_val,
              height=0.5,
              align='center',
              color=color,
              alpha=0.5,
              label=label+" Sum: %0.2f"%sv)

    ax.set_yticks(ytick_pos)
    ax.set_yticklabels(obj_keys, size='x-small')
    ax.set_ylim(-0.5, len(obj_keys) + 0.5)
    ax.set_xlabel('Objective value (# std)')
    ax.set_ylabel('Objectives')
    ax.yaxis.grid(True)

    ax2 = ax.twinx()
    ax2.set_yticks(ytick_pos)
    ax2.set_yticklabels(feat_val)
    ax2.set_ylim(-0.5, len(obj_keys) + 0.5)
    ax2.set_ylabel('Feat. values')
    ax.legend()
    
    figs[figname]['obj_keys'] = obj_keys
    fig.suptitle(reportname, fontsize=10)

    return fig


def diversity(checkpoint, evaluator, color='b', figs={}, reportname=''):
    '''plot the whole history, the hall of fame, and the best individual
    from a unumpyickled checkpoint
    '''

    param_names = evaluator.param_names
    n_params = len(param_names)

    hof = checkpoint['halloffame']
    fitness_cut_off = 2.*sum(hof[0].fitness.values)

    figname = 'Diversity ' + reportname
    fig = pt.make_figure(figname=figname,
                        orientation='page',
                        figs=figs,
                        fontsizes = (8,8))

    axs = pt.tiled_axs(frames=n_params, columns=6, d_out=5, fig=fig,
                top_margin=0.08, bottom_margin=0.03,
                left_margin=0.08, right_margin=0.03,
                hspace=0.3, wspace=0.8)

    all_params = get_params(checkpoint['history'].genealogy_history.values(),
                            fitness_cut_off=fitness_cut_off)
    hof_params = get_params(checkpoint['halloffame'],
                            fitness_cut_off=fitness_cut_off)

    best_params = checkpoint['halloffame'][0]

    param_count = len(param_names)
    for i, name in enumerate(param_names):

        ax = axs[i]

        p1 = numpy.array(all_params)[:, i].tolist()
        p2 = numpy.array(hof_params)[:, i].tolist()

        df = pd.DataFrame()
        df['val'] = p1+p2
        df['type'] = ['all']*len(p1)+['hof']*len(p2)
        df['param'] = [name]*(len(p1)+len(p2))

        sns.violinplot(x='param', y='val', hue='type', data=df, ax=ax, split=True,
               inner="quart", palette={"all": "gray", "hof": color})

        ax.axhline(y=best_params[i], color=color, label='best', linewidth=2)
        ax.set_ylabel('')

        if i > 0:
            ax.legend_.remove()
        else:
            ax.legend()


    for i, parameter in enumerate(evaluator.params):
        min_value = parameter.lower_bound
        max_value = parameter.upper_bound

        ax = axs[i]
        ax.set_ylim((min_value-0.02*min_value, max_value+0.02*max_value))

        name = param_names[i]
        label = name.replace('.', '\n')

        ax.set_title(label, fontsize=7)

        pt.adjust_spines(ax, ['left'], d_out=5)

        #ax.set_xticks([])
        #ax.axes.get_xaxis().set_visible(False)


    fig.suptitle(reportname, fontsize=10)

    return fig


def separate_params(
        axs,
        params,
        marker,
        color,
        markersize=40,
        fitness_cut_off=1e9,
        alpha=1):

    '''plot the individual parameter values'''
    observations_count = len(params)
    param_count = len(params[0])

    results = numpy.zeros((observations_count, param_count))
    good_fitness = 0
    for i, param in enumerate(params):
        if fitness_cut_off < sum(param.fitness.values):
            continue
        results[good_fitness] = param
        good_fitness += 1

    for i in range(param_count):
        ax = axs[i]
        #ax.set_clip_on(False)
        x = numpy.ones(good_fitness)
        y = results[:good_fitness, i]
        ax.scatter(x=x, y=y, s=float(markersize),
                    marker=marker, color=color,
                    clip_on=False,
                    alpha=alpha, edgecolors='none')

def get_params(
        params,
        fitness_cut_off=1e9
        ):

    '''plot the individual parameter values'''
    results = []
    for i, param in enumerate(params):
        if fitness_cut_off > sum(param.fitness.values):
            results.append(param)

    return results




def history(history):
    """Plot the history of the individuals"""

    import networkx

    plt.figure()

    graph = networkx.DiGraph(history.genealogy_tree)
    graph = graph.reverse()     # Make the grah top-down
    # colors = [\
    #        toolbox.evaluate(history.genealogy_history[i])[0] for i in graph]
    positions = networkx.graphviz_layout(graph, prog="dot")
    networkx.draw(graph, positions)
