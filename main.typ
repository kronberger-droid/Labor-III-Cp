#import "lab-report.typ": lab-report
#import "@preview/lilaq:0.2.0": (
    diagram,
  linspace,
  marks,
  path,
  place,
  plot,
  scatter,
  xaxis,
  yaxis,
)
#import "lib.typ": *
#import calc: ceil, pow

// create title page and initialize all of the formating
#show: lab-report.with(
  title: "Specific Heat of Solids",
  authors: ("Raul Wagner", "Martin Kronberger"),
  supervisor: "Juan Carlos Moreno Lopez",
  groupnumber: "301",
  date: datetime(day: 30, month: 4, year: 2025),
)

#set heading(numbering: "1.1")

#let resistance = 350 // in Ohm
#let mass = 4.995 // in g
#let molar_mass = 164.93 // g / mol
#let current = 25e-3 // A
#let measurement_current = 1e-3 //A

= Problem Statement

The specific heat, or specific heat capacity can be described as a material’s ability to absorb energy as a function of temperature.
Measuring specific heat provides experimental data essential for determining core thermodynamic properties such as entropy and the partition function.
It also offers insights into the elastic properties of solids.
Furthermore, the temperature dependence of specific heat reveals information about the system’s order, such as magnetic ordering or the onset of superconductivity.
In metals and alloys, in addition to the lattice-specific heat common to all solids, there is also an electromagnetic contribution, reflecting the behavior of conduction electrons and thus characterizing the metallic state.

= Experimental Setup

The experimental setup consists of a gas-fillable system containing a #mass g sample of Holmium. The lower part of the system containing the sample, is immersed in a Glass Dewar filled with liquid nitrogen. The sample is equipped with a Pt100 resistance thermometer used to measure its temperature and can be heated via a strain gauge ($350 plus.minus 1 space Omega $) which is attached to it.
To create nearly adiabatic conditions for the sample, the system can be evacuated to below 10⁻⁴ mbar using a rotary vane pump and a turbomolecular pump.
The wires used for measurement are additionally shielded as effectively as possible using copper plates and a copper shield surrounding the sample.
Furthermore, the sample is suspended on nylon threads to minimize additional heat conduction.
For cooling the sample, the system can be filled with helium through a valve, allowing thermal exchange between the nitrogen-cooled outer wall and the sample.
Separate power supplies are used for temperature measurement via the Pt100 resistance thermometer and for heating the sample via the strain gauge.
The Pt100 thermometer is powered with a current of 1 mA, and the resulting voltage drop is measured directly at the resistor using the four-point method.
The strain gauge is powered by a separate power source providing 25 mA.
Voltage measurements are carried out using the NI-6210 data acquisition module, which enables data evaluation in a dedicated software application.

#figure(
  image("assets/measurement_setup.png"),
  caption: [Setup: 1. exhaust gas, 2. rotary vane pump, 3. shut-off valve, 4. turbomolecular pump, 5. Penning vacuum gauge, 6. thermocouple vacuum gauge, 7. cryostat, 8. Dewar flask, 9. pressure gauge, 10. vent valve, 11. helium gas cylinder.]
)

== Cooling cycle before measurements <sec:cooling>

The sample was cooled down before the experiments following this procedure.
- Closing of the valve connecting the system and the rotary vane pump
- Switching off the turbomolecular pump and waiting for 15 minutes
- Switching off pressure gauges
- Opening the valve connecting the system to the helium gas supply and filling the system with 0,1 bar
- Closing the valve again and waiting until the desired temperature is reached
- Opening of the valve connecting the system with the rotary vane pump, switching on the thermocouple pressure gauge and waiting until the system has reached $10^(-2)$ mbar
- switching on the turbomolecular pump and the Penning vacuum gauge

As soon as the pressure is down to $10^(-4)$ mbar a measurement can be initiated. @Mueller2009

Since the speed of cooling down a material is proportional to the difference of the material and the environment temperature, it cools down rapidly in the beginning and slows down over time finally behaving asymptotic towards the surrounding temperature.

$
  (d T)/(d t) = -k(T_m -T_e)
$

== General Values of the System

#grid(
  columns: (auto, auto),
  align: (left, right),
  column-gutter: 10pt,
  gutter: 5pt,
  [Sample Mass:], [#mass $g / m$],
  [Sample Molar Mass:], [#molar_mass g],
  [Resistance of the DMS:], [$#resistance plus.minus 1 space Omega$],
  [Heating Current:], [$#(current * 1000) plus.minus 0.1$ mA],
  [Measurement Current (PT100)], [#(measurement_current * 1000) mA]
)

= Continuous Measurement

Using the Continuous Measurement Method to determine the specific heat capacity is a relatively simple method to get rough values for $C_p$ in a short time over a wide temperature range. In this method the sample is heated with a known amount of Energy determined by the heating current.
While heating the sample continuously the temperature was monitored through the NI-6210 data acquisition module via the LabVIEW Virtual Instrument program in an time interval of around half a second.
Since there is no thermal equilibrium between sample, heater and its wires during the measurement there are always unknown factors that can distort the measurement results, making this method inaccurate.

In this experiment the sample was cooled down to 80K according to the cooling procedure described in @sec:cooling.
Once the start temperature and a vacuum of $10^(-4)$ mbar was reached, the heating of the sample started with switching on the heater power supply using 25mA.
After reaching the final temperature of 300K the power-supply was turned off again and the experiment ended.

// parse data into a list of dicts:
// ( ( time: float, temperature: float, voltage: float), ... )
#let data = parse_measurements("data/spez.Warme.cont.lvm").slice(30, -400)

// approximation window width
#let window_width_low = 11;
#let window_width = 61;
#let window_width_high = 201;

To analyze the data acquired from the raw file provided by the measurement setup, it is parsed into data structure, consisting of a list of dictionaries:
#align(center)[
  ( ( time: float, temperature: float), ... )
]
After slicing of data points at the start and the end of the dataset which were not part of the experiment, but part of the resting periods between measurements, a window with a width of #window_width data-points is moved over the dataset.
At every location it uses a data-point close to the middle as time and temperature values corresponding to this point.
Afterwards $Delta T$ and $Delta t$ are approximated by assuming the curve inside of the current window to be a linear function through the boundary points.

#let calc_cont_cp_T(data, window_width) = {
  data.windows(window_width).map(window => {
    let delta_T = window.last().temperature - window.first().temperature
    let delta_t = window.last().time - window.first().time

    (
      time: window.at(ceil(window_width / 2)).time,
      temperature: window.at(ceil(window_width / 2)).temperature,
      cp: ((pow(current, 2) * resistance) / mass) * (delta_t / delta_T)
    )
  })
}

#let cp_T = calc_cont_cp_T(data, window_width)

#let cp_T_low = calc_cont_cp_T(data, window_width_low)

#let cp_T_high = calc_cont_cp_T(data, window_width_high)

This makes it possible to calculate the specific heat using the following formula for every position of the window, resulting in a near continuous curve.
$
  C_p (T) = 1 / m (d Q) / (d t) (d t) / (d T) = (I^2R) / m 1 / ((d T) / (d t))
$
where $m$ is the mass of the sample, $R$ is the resistance of the DMS and $I$ is the heating current.

// plot cp over T for the continuous measurement
#align(center)[
  #figure(
    diagram(
      width: 12cm,
      height: 8cm,
      title: [Continuous Measurement $C_p (T)$],
      xlabel: [Temperature $T$ in K],
      ylabel: [Specific heat $C_p$ in $J / (g K)$],
      plot(
        cp_T.map(x => x.temperature),
        cp_T.map(y => y.cp),
        label: [Window width: #window_width],
        mark: none,
        stroke: stroke(thickness: 1.5pt),
        color: blue
      ),
      plot(
        cp_T_low.map(x => x.temperature),
        cp_T_low.map(y => y.cp),
        label: [Window width: #window_width_low],
        mark: none,
        color: green.transparentize(60%),
        stroke: stroke(thickness: 0.8pt),
      ),
      plot(
        cp_T_high.map(x => x.temperature),
        cp_T_high.map(y => y.cp),
        label: [Window width: #window_width_high],
        mark: none,
        color: red.transparentize(60%),
        stroke: stroke(thickness: 0.8pt),
      ), 
    ),
    caption: [Specific heat capacity plotted over temperature acquired from continuous measurement method for different window widths.]
  )<fig:continuous>
]

The window width can be used to remove noise from the data, but results in loss of information particularly at locations of sudden changes in slope.
Like around the peak at approximately 130 K, where values get majorly distorted when choosing greater window widths.
Therefore we settled on a window width of #window_width which yields a relatively noise free dataset with still good precision.
It is also important to choose a odd window width to be able to map every window to the data point exactly in the middle of the window.
Otherwise this would have to shift the acquired data-points either one time-step to the right or left, resulting in unnecessary errors.

== Error analysis – Continuous heating
The Metering Time Period $t$ has 3 significant numbers with an error of +/- 0,0005s.

The measured temperature $T$ is measured by the resistance of the PT-100 resistance thermometer.
The error of the linearity between resistance and temperature is estimated with +/- 3%.
$
  delta T = plus.minus 0,03(T-27,781K)
$
/*Including the errors of all other values the error of our measured specific heat capacitance is the following:
$
  delta C_p = (0.05A dot 350Omega)/(#mass g) dot 1/((d T)/(d t)) dot 0.0005A + (0.625A)/(#mass g) dot 1/((d T)/(d t)) dot 1Omega + 
  (0.625A dot 350Omega)/(#mass g) dot 1/((d T)^2/(d t)) dot delta d T + 
  (0.625 A dot 350Omega)/(#mass g) dot 1/(d T) dot delta d t
$*/

Since the error increases with the temperature the maximum error is expected at the highest temperature.

Therefore the relative estimated error of $C_p$ can be $plus.minus 3,43%$ at 300K, taking the errors of all values into account. 


= Nerst Method

In order to get more accurate measurements, the Nerst Method is used.
In this method the sample is heated up in short bursts which are followed by short resting periods where no heat is supplied, allowing for a thermal equilibrium of the sample, the heater and its wires during the steps.
Measuring the temperature at this thermal equilibrium and taking thermal drifts of the sample, due to the not perfect adiabatic surroundings, into account, this method proves to be more accurate for measuring the specific heat capacity.

In this experiment the sample was again cooled down to 100K according to the cooling procedure described in @sec:cooling.
Once the start temperature and a vacuum of $10^(-4)$ mbar was reached, the heating of the sample started with switching on the heater power supply for 12 seconds using 25 mA.
The heat pulse was followed by a resting period of 8 seconds.
By repeating heating and resting pulse the heating of the sample continued in steps until it reached the end temperature of 170K. 
The temperature was monitored throughout the whole experiment using the NI-6210 data acquisition module and the LabVIEW Virtual Instrument program resulting in a dataset of temperature over time, with clearly observable heating and resting periods, marked by the voltage being high or low, as seen in the following plot:

// threshold voltage (max. approx 4 V min approx. 0 V):
#let threshold = 2.0

// prepare raw data by parsing the data from the file and thresholding voltage
// ( more information in the comments preceding the next plot)
#let step_data_test = threshold_data(
  parse_measurements("data/spez.Warme.step.lvm").slice(90, 150),
  threshold,
)

#let k_1 = -0.005
#let d_1 = 92.8

// plot a slice of the raw nerst method dataset for explanations
#align(center)[
  #figure(
    diagram(
      width: 10cm,
      height: 6cm,
      title: [Nerst method raw dataset slice T(t)],
      yaxis: (mirror: false),
      xlabel: [Time $t$ in s],
      ylabel: [Temperature $T$ in K],
      plot(
        step_data_test.map(x => x.time),
        step_data_test.map(y => y.temperature),
        mark: none,
      ),
      plot(
        linspace(45, 65, num: 100),
        linspace(45, 65, num: 100).map(y => (
          y * k_1 + d_1
        )),
        color: olive,
        stroke: stroke(thickness: 1pt, dash: "dashed"),
        mark: none,
      ),
      yaxis(
        position: right,
        label: [Heating voltage],
        ticks: (0, 1),
        plot(
          step_data_test.map(x => x.time),
          step_data_test.map(y => y.voltage),
          color: maroon,
          mark: none,
        )
      ),
    ),
    caption: [Raw slice of the data-set acquired using the Nerst method with thresholded voltage.]
  )
]

To be able to analyze and plot the given dataset, at first the voltage at every data point will be thresholded at #threshold V to make it easier to discern heating periods from non heating periods.
This thresholding actually results in a slight error since its possible that through applying the threshold data points can be attributed to the wrong set, leading to a uncertainty of a maximum of 1 time-step if both boundary points of a heating or resting cycle get thresholded to the wrong set.

// raw parsed data as list of dicts:
// ( time: float, temperature: float, voltage: float )
#let raw_step_data = parse_measurements("data/spez.Warme.step.lvm").slice(
  30,
  -30,
)

// data after thresholding the voltage:
// ( ( time: float, temperature: float, voltage: bool ), ... )
#let thresholded_step_data = threshold_data(raw_step_data, threshold)

Now the dataset is present as a list of dictionaries similar to the one for the continuous measurement, but this time including the voltage indicating the state of the heating element:
#align(center)[
  ( ( time: float, temperature: float, voltage: bool ), ... )
]
Now the data can be compartmentalized into heating and cooling periods, using the rising and falling flanks of the voltage.
Also for convenience the data sets get paired into tuples including the leading rest period with the followed heating period, which makes it easy to loop over the now much smaller data set of tuples:
#align(center)[
  ( ( resting-set, heating-set ), ... )
]

// split data at rising and falling flanks:
// ( ( points where voltage == 0 ), ( points where voltage: 1), ( point where voltage == 0 ), ... )
#let compartmenalized_step_data = compartmentalize_data(thresholded_step_data)

// combine every leading off-state with its following on state into a list of tuples:
// ( (points where voltage == 0 , points where voltage == 1 ), ... )
#let step_data = segmentize_data(compartmenalized_step_data)

Now using the following formulas the specific heat capacity $C_p$ and its corresponding temperature value $T_M$ can be calculated for every tuple.
$
  C_p (T_m) = (I_h^2 R t_h) / (m Delta T_x) quad "and" quad T_m = T_2 + k_1 t_h / 2 + (Delta T_x) / 2
$
Where $I_h$ is the heating current, $R$ is the resistance of the strain gauge, $t_h$ the duration of the heating period, $m$ the mass of the sample, $Delta T_x$ the change in temperature during the heating period, $T_2$ temperature value at the beginning of the heating period and $k_1$ the slope of the temperature over time during the resting period. 

// calculate a data-point for every on-off tuple using the given equations resulting in a list of dicts:
// ( ( time_m: float, temperature_m: float , cp: float), ... )
#let step_cp_T = step_data.map(segment => {
  // extract off period and slice of some of the front to remove rest heating
  let off_period = segment.first().slice(4)

  // calculate how long the heating was off
  let off_period_delta_t = (
    off_period.last().time - off_period.first().time
  )

  // use linear interpolation to calculate slope of off period
  let off_period_slope = linear_fit_slope(off_period)

  // extract onn period
  let on_period = segment.last()

  // calculate how long the heating was on
  let on_period_delta_t = (
    on_period.last().time - on_period.first().time
  )

  // calculate how much the temperature changed while heating was on
  let on_period_delta_T = (
    on_period.last().temperature - on_period.first().temperature
  )

  // calculate the temperature approx. in the middle of the heating period
  let middle_temp = (
    segment.last().first().temperature
      + off_period_slope * off_period_delta_t / 2
      + on_period_delta_T / 2
  )

  // write everything into the list of dicts using given formulas
  (
    time: segment.last().first().time + on_period_delta_t / 2,
    temperature: middle_temp,
    cp: (pow(current, 2) * resistance * on_period_delta_t)
      / (mass * on_period_delta_T),
  )
})

// slice data set of continuous measurement to compare it to nerst method
#let slice_cp_T = cp_T.slice(70, -200)

// shift the continous measurement until it closely matches the stepped one
#let shift = 0.01
#let shift_cp_T = slice_cp_T.map(point => {
  (
    temperature: point.temperature,
    cp: point.cp + shift,
  )
})

// plot cp over T
#align(center)[
  #figure(
    diagram(
      width: 12cm,
      height: 8cm,
      title: [Nerst Method $C_p (T_m)$],
      xlabel: [Temperature $T$ in K],
      ylabel: [Specific heat $C_p$ in $J / (g K)$],
      ylim: (0.18, 0.34),
      plot(
        shift_cp_T.map(x => x.temperature),
        shift_cp_T.map(y => y.cp),
        mark: none,
        color: olive,
        label: [Continuous method shifted by #shift $J / (g K)$],
      ),
      plot(
        step_cp_T.map(x => x.temperature),
        step_cp_T.map(y => {
          y.cp
        }),
        color: maroon,
        label: [Stepped measurement using Nerst method],
      ),
    ),
    caption: [Specific heat capacity over temperature using the Nerst method, compared to shifted data-set from continuous measurement as seen in @fig:continuous.]
  )
]

The plot shows clearly that, even tho the Nerst method is more accurate by including the intrinsic change in temperature which is present even without heating, it is much less precise.
This intrinsic temperature change results the values corresponding to the Nerst method being shifted up approximately $#shift J / (g K^2)$ compared to the continuous measurement, conducted before.

== Error estimation

The heat pulse length $t_H$ is approximately 12 seconds long.
With a metering time resolution of 0.5 seconds $t_H$ has an error of $plus.minus$0,5s.

The measured temperature $T$ is measured by the resistance of the PT-100 resistance thermometer.
The error of the linearity between resistance and temperature is estimated with $plus.minus$3%.
$
  delta T = plus.minus 0,03(T-27,781K)
$
Since the error increases with the temperature the maximum error is expected at the highest temperature.
Therefore the relative estimated error of $C_p$ can be $plus.minus$8,95% at 170K, taking the errors of all values into account.

== Magnetic Entropy and Neél Temperature

The magnetic contribution to the specific heat arises because, in magnetic solids, part of the internal energy is used to build up magnetization.
At the Neél Temperature $T_N$, the order of the aligned magnetic moments is completely destroyed.
At low temperatures, below the Neél Temperature, antiferromagnetic order can exist but vanishes at and above the Neél Temperature, turning the material typically into a paramagnet.

mag. Entropy:
$
  S = integral d T C_p / T
$
The entropy S can be determined from a $C_p/T$ vs. $T$ diagram by integrating the area under the curve.
To determine the magnetic entropy, the lattice contribution to the specific heat is approximately extrapolated linearly.
The area under the resulting linearized triangle is then calculated and corresponds to the entropy of the magnetic transition.

// define indices which approx. correspond to the magnetic peak
#let triangle_indices = (12, 21, 23)
#let triangle = triangle_indices.map(index => {
  (
    step_cp_T.at(index).temperature,
    step_cp_T.at(index).cp / step_cp_T.at(index).temperature,
  )
})

// create path to draw triangle in plot
#let triangle_segments = triangle.zip(
  triangle.slice(1, triangle.len()) + (triangle.at(0),),
)

// prepare for draw function
#let triangle_verts = (
  (triangle_segments.at(0).at(0),) + triangle_segments.map(seg => seg.at(1))
)

// define Nerst Temperature using the max. value of the magnetic peak
#let T_N = step_cp_T.at(21).temperature

// calculate the area of the previously defined triangle using the cross product
#let S_mag = triangle_area(triangle)

#let shift_T = 0.0001

// plot cp/T over T including the mag. entropy triangle and the point at nerst temperature
#align(center)[
  #figure(
    diagram(
      width: 12cm,
      height: 8cm,
      title: [Nerst Method $(C_p (T_m))/T_m$],
      xlabel: [Temperature $T$ in K],
      ylabel: [Specific heat over $T$ - $C_p / T$ in $J / (g K^2)$],
      ylim: (1.0e-3, 3.0e-3),
      legend: (position: top + right),
      plot(
        slice_cp_T.map(x => x.temperature),
        slice_cp_T.map(y => (y.cp / y.temperature) + shift_T),
        stroke: olive,
        mark: none,
        label: [Continuous measurement shifted by #shift_T $J / (g K^2)$],
      ),
      plot(
        step_cp_T.map(x => x.temperature),
        step_cp_T.map(y => (
          y.cp / y.temperature
        )),
        color: maroon,
        label: [Stepped measurement using Nerst method],
      ),
      path(
        ..triangle_verts,
        stroke: none,
        fill: eastern.transparentize(90%),
        closed: true,
        label: [Area determining magnetic entropy $S_(m)$],
      ),
      place(
        step_cp_T.at(triangle_indices.at(1)).temperature,
        step_cp_T.at(triangle_indices.at(1)).cp
          / step_cp_T.at(triangle_indices.at(1)).temperature,
        align: bottom + left,
        pad(.3em)[$C_p (T_N)$],
      ),
      place(48%, 50%, pad(.3em)[$S_(m)$]),
      scatter(
        triangle.map(i => i.first()),
        triangle.map(i => i.last()),
        color: eastern,
        size: 6pt,
      ),
    ),
    caption: [Specific heat over T plotted over T acquired from Nerst method compared to shifted data-set from continuous measurement. Including Nerst temperature and area determining magnetic entropy.]
  )
]

Therefore the Nerst temperature and the magnetic entropy yield with relatively high errors coming from the fact that the data-set is very noisy and the vertices are not strictly determined:
$
  T_N = #round(T_N, digits: 0) plus.minus 3" K"
  quad quad "and" quad quad
  S_("mag") = #(round(S_mag, digits: 3)*1000) plus.minus 3 ("mJ") / "K"
$
#pagebreak()

#outline(
  title: "Figures",
  target: figure.where(kind: image),
)
#bibliography("bibliography.bib")