---
title: Part 2 — All You Need to Know About Relational Hyper-Event Modeling
nav_order: 2
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

# All You Need to Know About Relational Hyper-Event Modeling

## 1. Core RHEMs

**Relational hyper event sequence:** realization of a **marked point process**.

**Sequence of undirected hyperevents** (Lerner, 2025):

$$E = \{e_i = (S_i, t_i),\ i=1,\dots,n\}$$

where:
- $S_i$ is a subset of the **vertex set** $V$, composed of the entities involved simultaneously in the $i$-th relational hyperevent, i.e. the **participants**;
- $t_i$: **time** of the $i$-th relational hyperevent.

*E.g. Les Misérables, whose chapters can be considered as hyperevents where actors co-appear:*

$$E = \{(\{\text{MY, NP, MB}\}, t_1),\ (\{\text{MY, ME, MB}\}, t_2),\ (\{\text{MY}\}, t_3),\ (\{\text{MY, ME, CL, GE, MC, MB}\}, t_4),\ \ldots\}$$

**Sequence of directed hyperevents** (Lerner, 2025; Boschi et al., 2025):

$$E = \{e_i = (S_i, R_i, t_i),\ i=1,\dots,n\}$$

where:
- $S_i$ is a subset of the **sender set** $V^S$, representing the group of senders of the $i$-th relational hyperevent;
- $R_i$ is a subset of the **receiver set** $V^R$, representing the group of receivers;
- $t_i$: **time** of the $i$-th relational hyperevent.

*E.g. Scientific innovation, where each paper publication is a hyperevent (a group of authors publishes a paper citing a set of papers):*

$$E = \{(\{2,3\}, \{A,B\}, t_1{=}1950),\ (\{1,2,4\}, \{A,C\}, t_2{=}1951),\ (\{4\}, \{A,C\}, t_3{=}1953)\}$$

> **Remark:** An undirected hyperevent can be represented as a directed hyperevent with an empty receiver set, i.e. $(S_i, t_i) = (S_i, R_i = \emptyset, t_i)$, also referred to as a **meeting**.

In both cases, $n$ is the **number of hyperevents** and $i$ the **index of a generic hyperevent**.

---

**Longitudinal network / multivariate counting process** (Perry & Zhu, 2013):

$$\{N_{SR}(t)\}_{S \subseteq V^S,\ R \subseteq V^R,\ t \geq 0}$$

where:
- $N_{SR}(t)$: number of interactions from $S$ to $R$ until time $t$;
- $V^S$: **set of senders**; $V^R$: **set of receivers**.

If senders/receivers belong to the same group, $V^S = V^R = V$ (**vertex set**), with $p = |V|$.

**Doob-Meyer decomposition theorem:**

$$N_{SR}(t) = \Lambda_{SR}(t) + M_{SR}(t) = \int_0^t \lambda_{SR}(u)\,\mathrm{d}u + M_{SR}(t)$$

where:
- $\Lambda_{SR}(t)$: **cumulative intensity process** (structural part);
- $M_{SR}(t)$: **Martingale residual process** (noisy part);
- $\lambda_{SR}(t)$: **hazard** (also: **rate**, **intensity**): instantaneous rate at which $(S,R)$ occurs at $t$.

---

**Relational hyper event model (RHEM)** (Lerner & Lomi, 2021, 2023):

$$\lambda_{SR}(t) = W_{SR}(t) \times \lambda_0(t) \times \exp\!\Big\{f\!\left(\boldsymbol{x}_{SR}(t), t\right)\Big\}$$

where:
- $W_{SR}(t)$: **risk indicator** — equal to 1 if $(S,R)$ is at risk at time $t$;
- $\lambda_0(t)$: **baseline hazard**;
- $f(\boldsymbol{x}_{SR}(t),t)$: **event log-rate contribution**. E.g. **linear contribution** $f(x_{SR}(t),t) = \beta \cdot x_{SR}(t)$:
  - $f(\cdot,t)$: **contribution function** (also: **effect**);
  - $\boldsymbol{x}_{SR}(t)$: edge-specific **covariates**;
  - $t$: **calendar time**.

$(S,R)$ satisfying $W_{SR}(t)=1$ compose the **risk set** $\mathcal{R}_t$.

---

### 1.1 Exogenous & endogenous drivers for undirected hyperevents

**Exogenous and endogenous drivers** are measured by means of **covariates** $x_S(t)$, which can be seen as aggregations or summaries of **attributes** $z_s(t)$ of nodes in subsets $S' \subseteq S$.

- **Exogenous covariates:** do not depend on the network of (previous) events. One can consider both the **main effect** of an exogenous feature and the effect of related **homophily**.

  *E.g. Gender. The covariate capturing the main effect can be computed as the average of the feature over nodes in $S$:*
  $$\text{avg.z}(S,t) = \sum_{s \in S} z_s(t) / |S|$$

  *Gender **dissimilarity** (a negative effect indicates **homophily**):*
  $$\text{diff.z}(S,t) = \sum_{\{s,s'\} \in \binom{S}{2}} \frac{|z_s(t) - z_{s'}(t)|}{\binom{|S|}{2}}$$

- **Endogenous covariates:** depend on the network of (previous) events and allow us to measure **relational mechanisms**.

  *E.g. **Subset repetition**: the average number of prior joint events over all subsets $S' \subseteq S$ of size $\rho$:*

  $$\begin{aligned}
  \text{hy.deg}(S', t) &= \sum_{e_i:\, t_i < t} \mathbf{1}(S' \subseteq S_i) \\
  \text{sub.rep}^{(\rho)}(S,t) &= \sum_{S' \in \binom{S}{\rho}} \frac{\text{hy.deg}(S', t)}{\binom{|S|}{\rho}}
  \end{aligned}$$

---

## 2. Non-linear RHEMs

### 2.1 Directed hyperevents

This section reviews examples of **directed hyperevents**, specifically referring to the paper publication hyperevent defined above. Data concerning directed hyperevents manifest in the form of $E = \{e_i = (S_i, R_i, t_i),\ i=1,\dots,n\}$.

### 2.2 Exogenous & endogenous drivers for directed hyperevents

**Endogenous covariates:** depend on the network of (previous) events and allow us to measure **relational mechanisms**.

*E.g. **Subset repetition** for directed hyperevents: the average number of prior joint events involving a subset of senders and a subset of receivers of a given size in the past event history.*

$$\begin{aligned}
\text{cite}^{\text{aut}-\text{pap}}(S, R, t) &= \sum_{t_i < t} \omega(t - t_i) \cdot \mathbf{1}_{\{S \subseteq S_i \cap R \subseteq R_i\}} \\
\text{subrep}^{(\rho,\ell)}(S, R, t) &= \sum_{(S', R') \in \binom{S}{\rho} \times \binom{R}{\ell}} \dfrac{\text{cite}^{\text{aut}-\text{pap}}(S', R', t)}{\binom{|S|}{\rho} \times \binom{|R|}{\ell}}
\end{aligned}$$

Subset repetition is a very general concept and can be extended well beyond the notion of citation.

### 2.3 Beyond Linearity

**What can break Linearity and Time-Homogeneity?**
- **Saturation:** individuals face cognitive overload and become unable to effectively process additional information.
- **Forgetfulness:** propensity to lose track of past details over time.
- **External changes:** shifts in the environment, technology, or social norms that alter how factors influence the system.

We adapt the notation from Part 1 to distinguish basis function indices when combining the two types of effects:

1. **Fixed linear effect (FLE):**
$$f^{\text{FLE}}\!\big(x_{SR}(t), t\big) = \beta \cdot x_{SR}(t)$$

2. **Time-varying effect (TVE):**
$$f^{\text{TVE}}\!\big(x_{SR}(t), t\big) = \alpha(t) \cdot x_{SR}(t), \qquad \alpha(t) = \sum_{l=1}^{L} \alpha_l\, a_l(t)$$

3. **Non-linear effect (NLE):**
$$f^{\text{NLE}}\!\big(x_{SR}(t)\big) = \sum_{q=1}^{Q} \beta_q\, b_q\!\big(x_{SR}(t)\big)$$

**Joint time-varying non-linear effect (TVNLE):** coefficients of basis evaluations of covariates are not fixed but vary over time.

$$f^{\text{TVNLE}}\!\left(x_{SR}(t), t\right) = \sum_{q=1}^Q \Bigg(\sum_{l=1}^L \alpha_{ql}\, a_l(t)\Bigg) b_q\!\left(x_{SR}(t)\right)$$

To improve **model fitting stability** when using these types of effects:

1. Apply a **monotonic transformation** to each covariate $x$:
$$\dot{x} = 1 - \exp\!\left(-\frac{x}{\epsilon}\right)$$
   - Choose $\epsilon$ so that the transformed covariate (evaluated for events) is as close as possible to a uniform distribution.
   - Apply the same transformation to both events and non-events.

2. Transform event times $t$ using the **empirical cumulative distribution function**:
$$\dot{t} = \mathrm{ecdf}_t(t)$$