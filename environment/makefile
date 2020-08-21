ifndef NAME
override NAME = virsorter
endif

main:
	sudo docker build -f dockerfile -t nakor/virus_extraction . \
	&& sudo docker push nakor/virus_extraction \
	&& rm -f ~/.singularity_images.cache/nakor-virus_extraction.img

tool:
	sudo docker build -f $(NAME)_dockerfile -t nakor/$(NAME) . \
	&& sudo docker push nakor/$(NAME) \
	&& rm -f ~/.singularity_images.cache/nakor-$(NAME).img
