#import calc: abs, odd, round


#let parse_measurements(path) = {
  let lines = read(path).split("\n").slice(2, -1)
  lines.map(line => {
    let parts = line
      .split("\t")
      .map(x => x.replace(",", ".").trim())
      .filter(x => x != "")
      .map(x => float(x))
    (
      time: parts.at(0),
      temperature: parts.at(1),
      voltage: parts.at(2),
    )
  })
}

#let threshold_data(data, threshold) = {
  data.map(point => (
    time: point.time,
    temperature: point.temperature,
    voltage: if point.voltage < threshold { 0 } else { 1 },
  ))
}

#let compartmentalize_data(data) = {
  let buff = data.first()
  let acc = (data.first(),)
  let data_comp = ()

  for p in data.slice(1) {
    if buff.voltage == p.voltage {
      acc.push(p)
      buff = p
    } else {
      data_comp.push(acc)
      acc = (p,)
    }
    buff = p
  }
  data_comp.push(acc)
  data_comp
}


#let segmentize_data(data_comp) = {
  data_comp
    .windows(2)
    .enumerate()
    .filter(it => not odd(it.at(0)))
    .map(it => it.at(1))
}

#let linear_fit_slope(data) = {
  let n = data.len()
  let sum_x = data.map(p => p.time).sum()
  let sum_y = data.map(p => p.temperature).sum()
  let sum_xy = data.map(p => p.time * p.temperature).sum()
  let sum_x2 = data.map(p => p.time * p.time).sum()
  (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x)
}

#let rolling_avg(data, width) = {
  step_cp_T
    .windows(width)
    .map(seg => {
      (
        time: seg.at(ceil(width / 2)).time,
        temperature: seg.at(ceil(width / 2)).temperature,
        cp: seg.map(entry => entry.cp).sum() / width,
      )
    })
}

#let triangle_area(p) = {
  let u = (
    p.at(1).at(0) - p.at(0).at(0),
    p.at(1).at(1) - p.at(0).at(1),
  )
  let v = (
    p.at(2).at(0) - p.at(0).at(0),
    p.at(2).at(1) - p.at(0).at(1),
  )
  0.5 * abs(u.at(0) * v.at(1) - u.at(1) * v.at(0))
}
