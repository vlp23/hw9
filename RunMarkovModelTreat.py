import ParameterClassesTreat as Q
import MarkovModelClassesTreat as MarkovCls
import SupportMarkovModel as SupportMarkov
import scr.SamplePathClasses as PathCls
import scr.FigureSupport as Figs


# create a cohort
cohort = MarkovCls.Cohort(
    id=0,
    therapy=Q.Therapies.MONO)


# simulate the cohort
simOutputsTreat= cohort.simulate()

# Question 5 & 6
# graph survival curve
PathCls.graph_sample_path(
    sample_path=simOutputsTreat.get_survival_curve(),
    title='Survival curve',
    x_label='Simulation time step',
    y_label='Number of alive patients'
    )

# graph histogram of survival times
Figs.graph_histogram(
    data=simOutputsTreat.get_survival_times(),
    title='Survival times of patients with HIV',
    x_label='Survival time (years)',
    y_label='Counts',
    bin_width=1
)



# print the outcomes of cohort with therapy
SupportMarkov.print_outcomes(simOutputsTreat, 'Combo therapy:')

#Question 7
numberstrokes = MarkovCls.PatientStateMonitor(parameters=Q.ParametersFixed(Q.Therapies.MONO))

print ('Count of the number of strokes:',(numberstrokes.developed_stroke()))
