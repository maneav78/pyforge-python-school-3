# How to submit homework
## Forking
Go to quantori's repository and fork it, creating a copy of repository in your own account.
![Forking 1](/images/submit_hw/forking.png)
![Forking 2](/images/submit_hw/forking-2.png)

## Cloning
Go to **forked** repository, press "code" button and copy the URL for it
![Clone URL](/images/submit_hw/cloning.png)

- Open the terminal in your project folder and clone the repository  
`git clone <your repo url>`  
`cd <repo name>`

- Make a new branch and switch to it  
`git switch --create <branch name>`

- Push the branch to remote repository  
`git push --set-upstream origin <branch name>`

## Adding reviewers
To add someone as reviewer, he should be collaborator in your repo.
In your **own** repo, open Settings - Collaborators and press "Add People"
![Collaborators 1](/images/submit_hw/collaborators.png)
![Collaborators 2](/images/submit_hw/collaborators-2.png)

## Editing
- Open project and make your own changes
- Stage the changes  
`git add .`
- Commit the changes
`git commit -m '<commit message>'`
- Push the commit to remote  
`git push`

## Pull request
New branch will appear on Github
![New Branch](/images/submit_hw/new_branch.png)

Open Pull Requests tab and Press "New Pull Request"
![Pull Request 1](/images/submit_hw/pull_request.png)
Change the target repo for pull request, you want to make PR to your **own** repo.
![Pull Request 2](/images/submit_hw/pull_request-2.png)
And change the incoming branch, from new branch to main.
![Pull Request 3](/images/submit_hw/pull_request-3.png)
Now create a pull request
![Pull Request 4](/images/submit_hw/pull_request-4.png)
Add reviewers (they must be collaborators in your repo) and create pull request.
![Pull Request 5](/images/submit_hw/pull_request-5.png)
