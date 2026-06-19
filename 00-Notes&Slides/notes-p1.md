---
title: Introduction to Relational Event Models
nav_order: 1
parent: Notes and Slides
math: mathjax
---

<style>
blockquote {
  background: #f0f7ff;
  border-left: 3px solid #93c5fd;
  border-radius: 6px;
  padding: 1rem 1.25rem;
  margin: 1.5rem 0;
  color: #111;
}
blockquote p { margin: 0.4rem 0; }
blockquote strong { color: #000; }
pre, code {
  background: #f8f8f8;
  border: 1px solid #e5e7eb;
  border-radius: 6px;
}
pre {
  padding: 1rem 1.25rem;
  overflow-x: auto;
}
code {
  padding: 2px 5px;
  font-size: 13px;
}
</style>

# Introduction to Relational Event Models

<div class="download-box">
  <strong>📥 Download files</strong>
  <ul>
    <li><a href="./Sunbelt-Daytona/Sunbelt26_Notes.pdf" download>Sunbelt26_Notes.pdf</a> — PDF version of the notes</li>
    <li><a href="" download>Sunbelt26_SlidesW1.pdf</a> — Slides</li>
  </ul>
</div>

<!-- "./Sunbelt-Daytona/Slides/01-Introduction-to-REMs/" -->

## 1. Core REMs

### Notation

**Relational event sequence:** realization of a **marked point process**.

$$E = \{e_i = (s_i, r_i, t_i),\ i=1,\dots,n\}$$

where:
- $s_i$: **sender** of the $i$-th relational event;
- $r_i$: **receiver** of the $i$-th relational event;
- $t_i$: **time** of the $i$-th relational event;
- $n$: **number of events**;
- $i$: **index of a generic event**.

---

### Longitudinal network / multivariate counting process
(Perry & Wolfe, 2013)

$$\{N_{sr}(t)\}_{s \in V^S,\ r \in V^R,\ t \geq 0}$$

where:
- $N_{sr}(t)$: number of interactions from $s$ to $r$ until time $t$;
- $V^S$: **set of senders** in the system;
- $V^R$: **set of receivers** in the system.

If senders/receivers belong to the same group, then $V^S = V^R = V$ (**vertex set**), with $p=\lvert V \rvert$.

**Doob-Meyer decomposition theorem:**

$$N_{sr}(t) = \Lambda_{sr}(t) + M_{sr}(t) = \int_{0}^t \lambda_{sr}(u)\,\mathrm{d}u + M_{sr}(t)$$

where:
- $\Lambda_{sr}(t)$: **cumulative intensity process** (structural part);
- $M_{sr}(t)$: **Martingale residual process** (noisy part);
- $\lambda_{sr}(t)$: **hazard** (also: **rate**, **intensity**): instantaneous rate at which $(s,r)$ occurs at $t$.

---

### Relational Event Model (REM)
(Bianchi et al., 2024)

$$\lambda_{sr}(t) = W_{sr}(t) \times \lambda_0(t) \times \exp\{f(\boldsymbol{x}_{sr}(t),t)\}$$

where:
- **Risk indicator** $W_{sr}(t)$ - equal to 1 if $(s,r)$ is at risk at time $t$;
- **Baseline hazard** $\lambda_0(t)$ (including global risk determinants);
- **Event log-rate contribution** $f(\boldsymbol{x_{sr}}(t),t)$. 
  E.g. **linear contribution**: $f(\boldsymbol{x_{sr}}(t),t) = \beta \cdot x_{sr}(t)$.
  - $f(\cdot,t)$: **contribution function** (also: **effect**);
  - $\boldsymbol{x}_{sr}(t)$: edge-specific **covariates**;
  - $t$: **calendar time**.

The dyads $(s,r)$ satisfying $W_{sr}(t)=1$ compose the **risk set** $\mathcal{R}_t$.

---

### Simulating Relational Event Data

**Gillespie algorithm for simulating REMs:**

$$\begin{aligned}
\Delta T_{sr}(t) &\sim \mathrm{Exponential}\!\left(\lambda_{sr}(t)\right), \\[6pt]
T = \min_{(s,r)\in\mathcal{R}_t} \Delta T_{sr}(t) &\sim \mathrm{Exponential}\!\left(\sum_{(s,r)\in\mathcal{R}_t} \lambda_{sr}(t)\right), \\[6pt]
D &\sim \mathrm{Multinomial}\!\left(\{p_{sr}(t)\}_{(s,r)\in\mathcal{R}_t}\right), \\[6pt]
p_{sr}(t) &= \frac{\lambda_{sr}(t)}{\sum_{(s,r)\in\mathcal{R}_t} \lambda_{sr}(t)}, \quad (s,r)\in\mathcal{R}_t.
\end{aligned}$$

where:
- $\Delta T_{sr}(t)$: **waiting time** until an event involving sender $s$ and receiver $r$ occurs, evaluated at time $t$;
- $T$: waiting time until the next event, evaluated at time $t$;
- $D$: identifies the dyad $(s,r)$ involved in the event occurring at time $t + T$.

> **Assumption:** $\lambda_{sr}(t)$ does not change between $t$ and $t + T$.

---

### Explanatory variables in REMs

**Exogenous & endogenous drivers:** measured by means of **covariates** $x_{sr}(t)$. Covariates can be seen as aggregations/summaries of information derived from **attributes** $z_{sr}(t)$ of nodes and edges.

- **Exogenous covariates:** do not depend on the network of (previous) events.

  *E.g. Work experience dissimilarity, obtained as the absolute value of the difference in work experience between the two nodes:*
  $$x_{sr}(t) = |WE_s - WE_r|$$

- **Endogenous covariates:** depend on the network of (previous) events and allow us to measure **relational mechanisms**.

  *E.g. Reciprocity.*

  - **Binary covariates:** check for the presence (in the previous event history) of **at least one event satisfying a condition** of interest.

    *Is there at least one reciprocal event in the past?*
    $$x_{sr}(t) = \mathbf{1}\!\big(\{(r,s,t^*) \mid t^* < t\} \neq \varnothing\big)$$

  - **Count covariates:** count the **number of events satisfying a condition** of interest.

    *How many reciprocal events occurred in the past?*
    $$x_{sr}(t) = \big|\{(r,s,t^*) \mid t^* < t\}\big|$$

  - **Time-aware covariates:** defined as a function of the **internal time** of the mechanism of interest (i.e., the difference between the current time and the time of the event satisfying the condition of interest).

    *How recent is the reciprocal event in the past?*
    $$z_{sr}(t) = \begin{cases} t - \tau_{sr}, & \tau_{sr} = \text{time of the last event } (r,s), \\ \infty, & \text{otherwise,} \end{cases}$$

    where $x_{sr}(t) = \exp\left(-z_{sr}(t)\right)$.

  Relational mechanisms can be **dyadic** or **triadic**, depending on the number of nodes involved. Some examples:
  - **Repetition** (dyadic)
  - **Transitive closure** (triadic)
  - **Cyclic closure** (triadic)

---

## 2. REM Inference

> **Assumption** (this section only): **linear contribution function**, i.e., 
$$f(\boldsymbol{x_{sr}}(t),t) = \boldsymbol{{\beta}} \cdot \boldsymbol{x_{sr}}(t)$$

The counting process is adapted to the **filtration**:

$$\mathbb{W} = \{\mathcal{W}_t\}_{t \in \mathbb{R}^+}$$

At time $t$, $\mathcal{W}_{t^-}$ incorporates both exogenous information and information from events occurred at $t^\ast < t$.

**Likelihood:** constructed as joint probability of observing data under the model.

Given a set of event times $t_1, t_2, \ldots, t_n$, which are realizations of a point process, the joint probability density of observing exactly $n$ event times in the interval $(0, t_n)$ is:

$$f(t_1, t_2, \ldots, t_n) = \prod_{i=1}^n \lambda(t_i \mid \mathcal{W}_{t_i^-}) \exp\left\{-\int_{t_{i-1}}^{t_i} \lambda(u \mid \mathcal{W}_{u^-})\,\mathrm{d}u\right\}$$

### 2.1 Full likelihood approaches

**Full likelihood** (Butts, 2008):

$$\mathcal{L}(\boldsymbol{\beta}) = \prod_{i=1}^n \underbrace{\sum_{(s,r) \in \mathcal{R}_{t_i}} \lambda_{sr}(t_i)\exp\!\left[-\sum_{(s,r) \in \mathcal{R}_{t_i}} \int_{t_{i-1}}^{t_i} \lambda_{sr}(u)\,\mathrm{d}u\right]}_{\text{Minimum Generalized Exponential Inter-arrival Time}} \underbrace{\dfrac{\lambda_{s_ir_i}(t_i)}{\sum_{(s,r) \in \mathcal{R}_{t_i}} \lambda_{sr}(t_i)}}_{\text{Multinomial Probability}}$$

**Piecewise Constant Rates (PCR) assumption** (Stadtfeld & Block, 2017): event rate is assumed to be constant within inter-arrival times, giving Exponential **waiting time** $\Delta T_i = T_i - T_{i-1}$.

$$\mathcal{L}(\boldsymbol{\beta}) \overset{\text{PCR}}{=} \prod_{i=1}^n \lambda_{s_ir_i}(t_i)\prod_{(s,r) \in \mathcal{R}_{t_i}} \exp\!\Bigg(-\lambda_{sr}(t_i)\left(t_i - t_{i-1}\right)\Bigg)$$

Letting $\Delta N_{s_ir_i} = N_{s_ir_i}(t_i) - N_{s_ir_i}(t_{i-1})$, this is proportional to a **Poisson Regression** likelihood.

> **Full log-likelihood in practice (PCR)**
>
> **Log-likelihood:**
> $$\ell(\boldsymbol{\beta}) = \sum_{i=1}^n \sum_{(s,r) \in \mathcal{R}_{t_i}} \Delta N_{sr}(t_i) \log\!\left(\lambda_{sr}(t_i;\boldsymbol{\beta})\right) - \lambda_{sr}(t_i;\boldsymbol{\beta})(t_i - t_{i-1})$$
>
> **Inference technique:** Poisson regression:
> $$\Delta N_{sr}(t_i)\mid \boldsymbol{x}_{sr}(t_i) \stackrel{\text{iid}}{\sim} \text{Poisson}\!\left(\mu_{sr}(t_i)\right)$$ 
> $$\log(\mu_{sr}(t_i)) = \boldsymbol{\beta} \cdot \boldsymbol{x}_{sr}(t_i) + \log(t_i - t_{i-1})$$
>
> **R code:**
> ```r
> glm(event ~ log_dist + offset(logtm),
>     family = poisson(link = "log"), data = poisson_dat)
> ```

### 2.2 Partial likelihood approaches

The **partial likelihood** (Cox, 1975) is the joint product of the multinomial terms:

$$\mathcal{L}^P(\boldsymbol{\beta}) = \prod_{i=1}^n\dfrac{\lambda_{s_ir_i}(t_i;\boldsymbol{\beta})}{\sum_{(s,r) \in \mathcal{R}_{t_i}} \lambda_{sr}(t_i;\boldsymbol{\beta})} = \prod_{i=1}^n \dfrac{\exp\{\boldsymbol{\beta}\cdot\boldsymbol{x}_{s_ir_i}(t_i)\}}{\sum_{(s,r) \in \mathcal{R}_{t_i}} \exp\{\boldsymbol{\beta}\cdot\boldsymbol{x}_{sr}(t_i)\}}$$

> **Partial log-likelihood in practice**
>
> **Log-likelihood:**
> $$\ell^P(\boldsymbol{\beta}) = \sum_{i=1}^n \left[\boldsymbol{\beta}\cdot\boldsymbol{x}_{s_ir_i}(t_i) - \log\!\left(\sum_{(s,r) \in \mathcal{R}_{t_i}} \exp\{\boldsymbol{\beta}\cdot\boldsymbol{x}_{sr}(t_i)\}\right)\right]$$
>
> **Inference technique:** CoxPH regression:
> $$\hat{\boldsymbol{\beta}} = \arg\max_{\boldsymbol{\beta}}\ell^P(\boldsymbol{\beta})$$
>
> **R code:**
> ```r
> rem(event ~ dist, data = cox_dat, method = "clogit", stratum = "stratum")
> ```

### 2.3 Sampled partial likelihood approaches

**Nested case-control dataset** (Borgan et al., 1995; Vu et al., 2015; Lerner & Lomi, 2020):

$$D_m = \{(e_i, \mathcal{SR}_{t_i}),\ \mathcal{SR}_{t_i} \subset \mathcal{R}_{t_i},\ |\mathcal{SR}_{t_i}|=m+1,\ i=1,\dots,n\}$$

where:
- $\mathcal{SR}_{t_i}$: **sampled risk set** at $t_i$, including $(s_i, r_i)$ and $m$ randomly sampled dyads in $\mathcal{R}_{t_i}$;
- $m$: **number of sampled non-events**.

**Relational event model with sampling:**

$$\lambda_{sr,\mathcal{SR}}(t) = \lambda_{sr}(t) \cdot \pi_t(\mathcal{SR}|(s,r)), \qquad \pi_t(\mathcal{SR}|(s,r)) = \frac{m}{|\mathcal{R}_t|-1} \cdot \mathbf{1}_{\{(s,r) \in \mathcal{SR},\, \mathcal{SR} \subset \mathcal{R}_t\}}$$

**Sampled partial likelihood:**

$$\mathcal{L}^{S}(\boldsymbol{\beta}) = \prod_{i=1}^n \dfrac{\exp\{\boldsymbol{\beta}\cdot\boldsymbol{x}_{s_ir_i}(t_i)\}}{\sum_{(s,r) \in \mathcal{SR}_{t_i}} \exp\{\boldsymbol{\beta}\cdot\boldsymbol{x}_{sr}(t_i)\}}$$

which coincides with the likelihood of a **conditional logistic regression** of a stratum of size $m+1$.

> **Sampled partial log-likelihood in practice**
>
> **Log-likelihood:**
> $$\ell^S(\boldsymbol{\beta}) = \sum_{i=1}^n \Bigg[\boldsymbol{\beta}\cdot\boldsymbol{x}_{s_ir_i}(t_i) - \log\!\Bigg(\sum_{(s,r) \in \mathcal{SR}_{t_i}} \exp\{\boldsymbol{\beta}\cdot\boldsymbol{x}_{sr}(t_i)\}\Bigg)\Bigg]$$
>
> **Inference technique:** Conditional logistic regression:
> $$\Delta N_{sr}(t_i)\,\Big|\sum_{s'r' \in \mathcal{SR}_i} \Delta N_{s'r'}(t_i) = 1 \sim \text{Categorical}\!\left(\pi_{sr}(t_i)\right)$$ 
>
>where:
>$$\pi_{sr}(t_i) = \dfrac{\exp\{\boldsymbol{\beta}\cdot\boldsymbol{x}_{sr}(t_i)\}}{\sum_{(s',r') \in \mathcal{SR}_{t_i}} \exp\{\boldsymbol{\beta}\cdot\boldsymbol{x}_{s'r'}(t_i)\}}$$
>
> **R code:**
> ```r
> rem(event ~ dist, data = sampled_dat, method = "clogit", stratum = "stratum")
> ```

### 2.4 Case-1-control partial likelihood approaches

Consider the sampled partial likelihood in the particular case of $m=1$ (Boschi et al., 2025; Filippi et al., 2024):

$$\mathcal{L}^{S}(\boldsymbol{\beta}) \stackrel{m=1}{=} \prod_{i=1}^n \operatorname{logistic}\!\big(\boldsymbol{\beta} \cdot \Delta\boldsymbol{x}_i\big)$$

where:
- $(s_i^\ast, r_i^\ast)$: the **sampled dyad** at time $t_i$, drawn from the risk set $\mathcal{R}_{t_i}$;
- $\Delta\boldsymbol{x}_i = \boldsymbol{x}_{s_ir_i}(t_i) - \boldsymbol{x}_{s_i^\ast r_i^\ast}(t_i)$: the difference between covariates for the event and the sampled non-event.

This coincides with the likelihood of a **degenerate logistic regression**: without intercept term, with constant response equal to 1, and with covariates equal to the difference $\Delta\boldsymbol{x}_i$.

> **Case-1-control likelihood in practice**
>
> **Log-likelihood:**
> $$\ell^S(\boldsymbol{\beta}) \stackrel{m=1}{=} \sum_{i=1}^n \Bigg[\boldsymbol{\beta}\cdot\boldsymbol{x}_{s_ir_i}(t_i) - \log\!\Bigg(\sum_{(s,r) \in \{(s_i,r_i),(s_i^\ast r_i^\ast)\}} \exp\{\boldsymbol{\beta}\cdot\boldsymbol{x}_{sr}(t_i)\}\Bigg)\Bigg]$$
>
> **Inference technique:** Degenerate logistic regression:
> $$Y_i|\Delta\boldsymbol{x}_i \stackrel{\text{iid}}{\sim} \text{Bernoulli}(\pi_i),\quad \text{logit}(\pi_i) = \boldsymbol{\beta}\cdot\Delta\boldsymbol{x}_i,\quad y_i = 1,\ i=1,\ldots,n;\quad \text{no intercept}$$
>
> **R code:**
> ```r
> rem(~ dist, data = dat_gam, method = "gam")
> ```

---

## 3. Non-Linear REMs

**Additive relational event model** (example): consider the REM with

$$f\!\left(\boldsymbol{x}_{sr}(t), \boldsymbol{v}_{sr}(t), t\right) = \underbrace{\beta \cdot x_{sr}^{(1)}(t)}_{\text{fixed linear}} + \underbrace{\beta(t) \cdot x_{sr}^{(2)}(t)}_{\text{time-varying}} + \underbrace{f^{(3)}(x_{sr}^{(3)}(t))}_{\text{non-linear}} + \underbrace{\boldsymbol{\gamma} \cdot \boldsymbol{v}_{sr}(t)}_{\text{random}}$$

including four types of effects:
1. **Fixed linear effect (FLE):** $\beta$ is time-homogeneous and is not a function of $x_{sr}^{(1)}(t)$;
2. **Time-varying effect (TVE):** $\beta(t)$ is a non-linear function of time;
3. **Non-linear effect (NLE):** $f^{(3)}(x_{sr}^{(3)}(t))$ is a non-linear function of $x_{sr}^{(3)}(t)$;
4. **Random effect (RE):** $\boldsymbol{\gamma} = (\gamma_s)_{s \in V}$, $\gamma_s$ is a Gaussian random variable.

This model can be fitted using a degenerate logistic regression.

### Time-varying effect

In its simplest form, a **smooth function of time** can be represented as a linear combination of basis functions of **calendar time**:

$$\beta(t) \approx \sum_{j=1}^q \theta_j b_j(t)$$

where:
- $q$: **basis dimension**;
- $\{b_j\}_j$: set of **basis functions**;
- $\theta_j$: **basis coefficients**.

> **TVE in practice** (Boschi et al., 2025)
>
> **Contribution to the log-event rate:**
> $$f\!\left(\boldsymbol{x}_{sr}(t),\boldsymbol{v}_{sr}(t),t\right) = \cdots + \underbrace{\sum_{j=1}^q \theta_j b_j(t)}_{\text{TVE}} \times x_{sr}^{(2)}(t) + \cdots$$
>
> **R code:**
> ```r
> rem(~ tv(distPT), time = "time", data = dat_gam_tv, method = "gam")
> ```

### Non-linear effect

A **smooth function of a covariate** can be represented as a linear combination of basis functions of the covariate itself:

$$f(x_{sr}(t)) \approx \sum_{j=1}^q \theta_j b_j(x_{sr}(t))$$

> **NLE in practice** (Filippi et al., 2024)
>
> **Contribution to the log-event rate:**
> $$f\!\left(\boldsymbol{x}_{sr}(t),\boldsymbol{v}_{sr}(t),t\right) = \cdots + \underbrace{\sum_{j=1}^q \theta_j b_j(x_{sr}^{(3)}(t))}_{\text{NLE}} + \cdots$$
>
> **R code:**
> ```r
> rem(~ nl(dist), data = dat_gam_nl, method = "gam")
> ```

### Random effect

Can be efficiently estimated as **0-dimensional splines**, with basis functions taking the value 1 when the level of the random factor is present and 0 otherwise:

$$\sum_{s' \in V} \gamma_{s'}\mathbb{1}\{s' = s\}, \qquad \boldsymbol{\gamma} = (\gamma_1,\ldots,\gamma_{|V|})$$

> **RE in practice** (Boschi et al., 2025)
>
> **Contribution to the log-event rate:**
> $$f\!\left(\boldsymbol{x}_{sr}(t),\boldsymbol{v}_{sr}(t),t\right) = \cdots + \underbrace{\sum_{s' \in V} \gamma_{s'}\mathbb{1}\{s = s'\}}_{\text{RE}} + \cdots$$
>
> **R code:**
> ```r
> rem(~ re(sender), data = dat_gam_r, method = "gam")
> ```